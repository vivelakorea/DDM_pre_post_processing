import os
import re
import numpy as np
import json
import argparse
import subprocess
from functools import partial
from multiprocessing import Pool, cpu_count
import sys
import zipfile

# ==============================================================================
# ==                    Common Utilities & Path Setup                         ==
# ==============================================================================
CWD = os.path.abspath(os.path.dirname(__file__))
ABAQUS_DIR = os.path.join(CWD, "ABAQUS")
RESTART_DIR = os.path.join(CWD, "restart")
PARADIS_VTK_DIR = os.path.join(CWD, "vtk_paradis")
ZIP_OUTPUT_DIR = os.path.join(CWD, "zip_archives")
ABAQUS_WORKER_SCRIPT_NAME = "temp_abaqus_worker.py"

def execute_command(command):
    """Executes a shell command and checks for success, printing all output."""
    print(f"\n[EXEC] {command}", flush=True)
    result = subprocess.run(command, shell=True, capture_output=True, text=True, encoding='utf-8', errors='ignore')
    stdout = result.stdout.strip() if result.stdout else ""
    stderr = result.stderr.strip() if result.stderr else ""
    if stdout: print(f"--- stdout ---\n{stdout}", flush=True)
    if stderr: print(f"--- stderr ---\n{stderr}", flush=True)
    if result.returncode != 0 or "error" in stderr.lower():
        print(f"[ERROR] Command failed with exit code {result.returncode} or produced an error.")
        if 'abaqus' in command.lower(): print("[HINT] Ensure 'abaqus' is in your system's PATH.")
        return False
    print("[SUCCESS] Command executed successfully.")
    return True

def create_zip_archive(source_dirs, zip_filename):
    """Creates a zip archive from all relevant files in the source directories."""
    zip_path = os.path.join(ZIP_OUTPUT_DIR, zip_filename)
    print(f"\n--- Creating Zip Archive: {zip_filename} ---")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for source_dir in source_dirs:
            if not os.path.isdir(source_dir):
                print(f"[WARN] Source directory not found, skipping: {source_dir}")
                continue
            for root, _, files in os.walk(source_dir):
                for file in files:
                    if file.endswith(('.vtu', '.vtp', '.pvd')):
                        file_path = os.path.join(root, file)
                        archive_name = os.path.relpath(file_path, CWD)
                        print(f"  Adding: {archive_name}")
                        zipf.write(file_path, arcname=archive_name)
    print(f"[SUCCESS] Archive created: {zip_path}")

# ==============================================================================
# ==                      Abaqus Mode Functions                             ==
# ==============================================================================

def create_abaqus_worker_script():
    """Generates the Python script that will be executed by 'abaqus python'."""
    script_content = """
import sys
import os
import json
import argparse
import numpy as np
try:
    from odbAccess import openOdb
    from abaqusConstants import *
except ImportError:
    print("Failed to import Abaqus modules. This script must be run with 'abaqus python'.")
    sys.exit(1)

# --- Content from utilities.py ---
class ReadableOdb:
    def __init__(self, odb, readOnly=True): self._odb = openOdb(odb, readOnly)
    def getFrames(self, stepName): return self._odb.steps[stepName].frames
    def getFrame(self, stepName, frameNum): return self.getFrames(stepName)[frameNum]
    def getInstance(self, instanceName): return self._odb.rootAssembly.instances[instanceName]
    def getNodes(self, instanceName): return self.getInstance(instanceName).nodes
    def getElements(self, instanceName): return self.getInstance(instanceName).elements
    def getFieldOutputs(self, stepName, frameIdx): return self._odb.steps[stepName].frames[frameIdx].fieldOutputs.items()
    @property
    def getStepsKeys(self): return self._odb.steps.keys()
    @property
    def getInstancesKeys(self): return self._odb.rootAssembly.instances.keys()

# --- Content from odb2vtk.py ---
def ABAQUS_VTK_CELL_MAP(abaqusElementType):
    type_map = {'C3D4': 10, 'C3D6': 13, 'C3D8': 12, 'C3D8R': 12, 'C3D10': 24,
                'C3D15': 26, 'C3D20': 25, 'S3': 5, 'S4': 9, 'S8': 23,
                'S9': 28, 'B31': 3, 'R3D3': 5, 'R3D4': 9}
    for key, value in type_map.items():
        if key in abaqusElementType: return value
    return 1 # Fallback to VTK_VERTEX

def ABAQUS_VTK_FIELDOUPUTS_MAP(fldOutput):
    abaqusDataType = fldOutput.type
    abaqusComponentLabels = fldOutput.componentLabels
    abaqusPosition = fldOutput.locations[0].position
    vtkType = 'Tensors'
    if abaqusDataType == SCALAR:
        vtkType = 'Scalars'
        if len(abaqusComponentLabels) == 0: abaqusComponentLabels = ('0',)
    elif abaqusDataType == VECTOR:
        vtkType = 'Vectors'
    return (vtkType, abaqusComponentLabels, abaqusPosition)

class ODB2VTK:
    def __init__(self, fileFullName, outputDir, pvd_suffix=''):
        self.fileFullName = fileFullName
        self.outputDir = outputDir
        self.pvd_suffix = pvd_suffix
        self.odb = ReadableOdb(self.fileFullName)
        self._nodes_map, self._elements_map, self._instance_names, self._step_frame_map = {}, {}, [], {}
        self._nodesNum, self._cellsNum = 0, 0

    def ExtractHeader(self):
        dictJson = {"instances": list(self.odb.getInstancesKeys), "steps": []}
        for stepName in self.odb.getStepsKeys:
            frames = self.odb.getFrames(stepName)
            frame_names = [f"{stepName}-frame-{i:03d}" for i in range(len(frames))]
            dictJson['steps'].append([stepName, frame_names])
        json_path = self.fileFullName.replace(".odb", ".json")
        with open(json_path, 'w') as fp:
            json.dump(dictJson, fp, indent=4)
        print(f"Header written to {os.path.basename(json_path)}")
    
    def ReadArgs(self, instanceNames, stepsFramesDict): self._instance_names = instanceNames; self._step_frame_map = stepsFramesDict
    def ConstructMap(self):
        indexNode, indexElement = 0, 0
        for instanceName in self._instance_names:
            self._nodes_map[instanceName], self._elements_map[instanceName] = {}, {}
            nodes, elements = self.odb.getNodes(instanceName), self.odb.getElements(instanceName)
            self._nodesNum += len(nodes); self._cellsNum += len(elements)
            for node in nodes: self._nodes_map[instanceName][node.label] = indexNode; indexNode += 1
            for elem in elements: self._elements_map[instanceName][elem.label] = indexElement; indexElement += 1
            
    def WriteFieldOutputData(self, fldOutput, pointdata_map, celldata_map):
        if not fldOutput.values: return ('', '')
        vtkData = ABAQUS_VTK_FIELDOUPUTS_MAP(fldOutput)
        position = vtkData[2]
        if position == NODAL: return (self.WriteDataArray(fldOutput, vtkData, fldOutput.name, "PointData", pointdata_map), '')
        elif position in [INTEGRATION_POINT, CENTROID]: return ('', self.WriteDataArray(fldOutput, vtkData, fldOutput.name, "CellData", celldata_map))
        return ('', '')
        
    def WriteDataArray(self, fldOutput, vtkData, description, dataType, data_map):
        vtkType, labels, position = vtkData
        num_components = len(labels)
        buffer = f'<DataArray type="Float32" Name="{description}" NumberOfComponents="{num_components}" format="ascii">\\n'
        size = self._nodesNum if dataType == "PointData" else self._cellsNum
        dataArray = np.zeros((size, num_components))
        
        for instanceName in self._instance_names:
            try:
                subset = fldOutput.getSubset(region=self.odb.getInstance(instanceName), position=position)
                if not subset.values: continue
            except: continue
            
            data_values = np.array([v.data for v in subset.values])
            if num_components == 1 and data_values.ndim == 1:
                data_values = data_values.reshape(-1, 1)

            if dataType == "PointData":
                indices = [self._nodes_map[instanceName][v.nodeLabel] for v in subset.values]
                dataArray[indices] = data_values
            elif dataType == "CellData":
                elem_data = {v.elementLabel: data_values[i] for i, v in enumerate(subset.values)}
                valid_elements = [e for e in self.odb.getElements(instanceName) if e.label in elem_data]
                if valid_elements:
                    indices = [self._elements_map[instanceName][e.label] for e in valid_elements]
                    data_rows = [elem_data[e.label] for e in valid_elements]
                    if indices: dataArray[indices] = np.array(data_rows)

        buffer += "\\n".join([" ".join(map(str, row)) for row in dataArray])
        buffer += '\\n</DataArray>\\n'
        data_map.setdefault(vtkType, []).append(description)
        return buffer

    def WriteVTUFile(self, args):
        stepName, frameIdx = args
        os.makedirs(self.outputDir, exist_ok=True)
        vtu_path = os.path.join(self.outputDir, f"{stepName}_{frameIdx}.vtu")
        with open(vtu_path, 'w') as f:
            f.write('<VTKFile type="UnstructuredGrid" version="1.0">\\n'
                    f'<UnstructuredGrid>\\n<Piece NumberOfPoints="{self._nodesNum}" NumberOfCells="{self._cellsNum}">\\n')
            f.write('<Points>\\n<DataArray type="Float32" NumberOfComponents="3" format="ascii">\\n')
            for instanceName in self._instance_names:
                for node in self.odb.getNodes(instanceName):
                    f.write(f"{node.coordinates[0]} {node.coordinates[1]} {node.coordinates[2]}\\n")
            f.write('</DataArray>\\n</Points>\\n')
            connectivity, offsets, types, offset = [], [], [], 0
            for instanceName in self._instance_names:
                for cell in self.odb.getElements(instanceName):
                    connectivity.append(" ".join(str(self._nodes_map[instanceName][n]) for n in cell.connectivity))
                    offset += len(cell.connectivity)
                    offsets.append(str(offset))
                    types.append(str(ABAQUS_VTK_CELL_MAP(cell.type)))
            f.write('<Cells>\\n'
                    '<DataArray type="Int32" Name="connectivity" format="ascii">\\n' + "\\n".join(connectivity) + '\\n</DataArray>\\n'
                    '<DataArray type="Int32" Name="offsets" format="ascii">\\n' + "\\n".join(offsets) + '\\n</DataArray>\\n'
                    '<DataArray type="UInt8" Name="types" format="ascii">\\n' + "\\n".join(types) + '\\n</DataArray>\\n'
                    '</Cells>\\n')
            pointdata_map, celldata_map = {}, {}
            p_bufs, c_bufs = [], []
            for fldName, fldOutput in self.odb.getFieldOutputs(stepName, frameIdx):
                p_buf, c_buf = self.WriteFieldOutputData(fldOutput, pointdata_map, celldata_map)
                if p_buf: p_bufs.append(p_buf)
                if c_buf: c_bufs.append(c_buf)
            if p_bufs: f.write('<PointData>\\n' + "".join(p_bufs) + '</PointData>\\n')
            if c_bufs: f.write('<CellData>\\n' + "".join(c_bufs) + '</CellData>\\n')
            f.write('</Piece>\\n</UnstructuredGrid>\\n</VTKFile>')
        print(f"VTU file written: {os.path.basename(vtu_path)}")

    def WritePVDFile(self):
        odb_base_name = os.path.splitext(os.path.basename(self.fileFullName))[0]
        pvd_filename = f"{odb_base_name}{self.pvd_suffix}.pvd"
        pvd_path = os.path.join(self.outputDir, pvd_filename)
        with open(pvd_path, 'w') as f:
            f.write('<VTKFile type="Collection" version="1.0">\\n<Collection>\\n')
            for stepName, frameList in self._step_frame_map.items():
                for frameIdx in frameList:
                    frame = self.odb.getFrame(stepName, frameIdx)
                    vtu_filename = f"{stepName}_{frameIdx}.vtu"
                    f.write(f'    <DataSet timestep="{frame.frameValue}" file="{vtu_filename}"/>\\n')
            f.write('  </Collection>\\n</VTKFile>')
        print(f"PVD file written: {os.path.basename(pvd_path)}")

def main_abaqus_worker():
    parser = argparse.ArgumentParser(description="Internal Abaqus worker script.")
    parser.add_argument("worker_mode", choices=['header', 'odb2vtk', 'pvd'])
    parser.add_argument("--odbFile", required=True)
    parser.add_argument("--outputDir", help="The absolute path to the output directory.")
    parser.add_argument("--instance", type=str)
    parser.add_argument("--step", type=str)
    parser.add_argument("--pvd-suffix", default='')
    args = parser.parse_args(sys.argv[1:])

    outputDir = args.outputDir if args.outputDir else os.path.dirname(args.odbFile)
    odb2vtk = ODB2VTK(args.odbFile, outputDir, args.pvd_suffix)
    if args.worker_mode == 'header':
        odb2vtk.ExtractHeader()
    elif args.worker_mode == 'pvd':
        step_name, frame_indices_str = args.step.split(':')
        frame_indices = [int(i) for i in frame_indices_str.split(',')]
        odb2vtk.ReadArgs([args.instance], {step_name: frame_indices})
        odb2vtk.WritePVDFile()
    elif args.worker_mode == 'odb2vtk':
        step_name, frame_idx_str = args.step.split(':')
        frame_idx = int(frame_idx_str)
        odb2vtk.ReadArgs([args.instance], {step_name: [frame_idx]})
        odb2vtk.ConstructMap()
        odb2vtk.WriteVTUFile([step_name, frame_idx])

if __name__ == "__main__":
    main_abaqus_worker()
"""
    worker_script_path = os.path.join(CWD, ABAQUS_WORKER_SCRIPT_NAME)
    with open(worker_script_path, "w") as f:
        f.write(script_content)
    print(f"Generated Abaqus worker script: {ABAQUS_WORKER_SCRIPT_NAME}")
    return worker_script_path

# A wrapper function for os.system that can be pickled by multiprocessing.Pool
def run_command_for_pool(command):
    """A simple wrapper for os.system to be used with multiprocessing.Pool"""
    exit_code = os.system(command)
    if exit_code != 0:
        # Raise an exception to make the pool aware of the failure.
        raise RuntimeError(f"Worker command '{command}' failed with exit code {exit_code}")
    return exit_code

def abaqus_main(args, worker_script_path):
    print("--- Running Abaqus to VTK Conversion ---")
    odb_file_path = os.path.abspath(os.path.join(ABAQUS_DIR, args.odbFile))
    json_file_path = odb_file_path.replace(".odb", ".json")

    output_dir_name = os.path.splitext(os.path.basename(odb_file_path))[0]
    abaqus_output_dir = os.path.join(ABAQUS_DIR, output_dir_name)
    os.makedirs(abaqus_output_dir, exist_ok=True)

    if not os.path.exists(odb_file_path):
        print(f"[FATAL] ODB file not found: {odb_file_path}"); return None, None

    if not os.path.exists(json_file_path):
        print(f"'{os.path.basename(json_file_path)}' not found. Generating header...")
        cmd = f'abaqus python "{worker_script_path}" header --odbFile "{odb_file_path}"'
        if not execute_command(cmd): return None, None
        if not os.path.exists(json_file_path):
            print(f"[FATAL] JSON file was not created. Aborting."); return None, None
    else:
        print(f"'{os.path.basename(json_file_path)}' already exists. Skipping header generation.")

    with open(json_file_path, "r") as f: data = json.load(f)
    try:
        instance_name, frames = data["instances"][0], data["steps"][0][1]
        total_frame_count = len(frames)
        step_name = data['steps'][0][0]

        if args.n_frames and args.n_frames < total_frame_count:
            selected_indices = np.linspace(0, total_frame_count - 1, args.n_frames, dtype=int)
        else:
            selected_indices = list(range(total_frame_count))

        frame_count = len(selected_indices)
        print(f"Instance: {instance_name}, Selected Frames: {frame_count}")
    except (KeyError, IndexError):
        print(f"[ERROR] Invalid JSON format in '{os.path.basename(json_file_path)}'."); return None, None

    commands = []
    for i in selected_indices:
        step_frame = f"{step_name}:{i}"
        cmd = (f'abaqus python "{worker_script_path}" odb2vtk --outputDir "{abaqus_output_dir}" '
               f'--odbFile "{odb_file_path}" --instance "{instance_name}" --step "{step_frame}"')
        commands.append(cmd)

    print(f"\nProcessing {frame_count} frames in parallel using {cpu_count()} cores...")
    try:
        with Pool(processes=cpu_count()) as pool:
            for result in pool.imap_unordered(run_command_for_pool, commands):
                pass
    except Exception as e:
        print(f"\n[FATAL] A parallel worker failed. Aborting Abaqus processing.")
        print(f"Error details: {e}")
        return None, None

    # cmd_pvd = (f'abaqus python "{worker_script_path}" pvd --outputDir "{abaqus_output_dir}" '
    #            f'--odbFile "{odb_file_path}" --instance "{instance_name}" '
    #            f'--step "{step_name}:{','.join(map(str, selected_indices))}")

    # if not execute_command(cmd_pvd): return None, None

    print("\n--- Abaqus to VTK Conversion Finished ---\n")
    return frame_count, abaqus_output_dir

# ==============================================================================
# ==                      ParaDiS Mode Functions                            ==
# ==============================================================================
def paradis_main(args, frame_count_from_abaqus=None):
    print("--- Running ParaDiS to VTK Conversion ---")
    num_frames = args.n_frames if args.n_frames else (
        frame_count_from_abaqus if frame_count_from_abaqus is not None else args.paradis_frames)

    start_idx = 2
    try:
        restart_files = [f for f in os.listdir(RESTART_DIR) if f.startswith('rs') and f.endswith('.data')]
        if not restart_files:
            print(f"[ERROR] No ParaDiS restart files found in '{RESTART_DIR}'.")
            return None
        last_file_num = max([int(re.findall(r'\d+', f)[0]) for f in restart_files])
        end_idx = last_file_num
    except (FileNotFoundError, IndexError, ValueError):
        print(f"[ERROR] Could not determine the last frame number in '{RESTART_DIR}'.")
        return None
    print(f"ParaDiS frame range: {start_idx} to {end_idx}. Sampling {num_frames} frames.")
    indices = np.linspace(start_idx, end_idx, num_frames, dtype=int)

    process_func = partial(process_single_dislocation, args=args)
    if args.no_parallel:
        results = list(map(process_func, indices))
    else:
        with Pool(processes=cpu_count()) as pool:
            results = pool.map(process_func, indices)
    valid_results = [r for r in results if r is not None]
    if not valid_results:
        print("[ERROR] No valid ParaDiS VTP files were generated.")
        return None
    pvd_path = os.path.join(PARADIS_VTK_DIR, "paradis_series.pvd")
    write_paradis_pvd(pvd_path, valid_results)
    print("\n--- ParaDiS to VTK Conversion Finished ---\n")
    return PARADIS_VTK_DIR


def parse_params(param_file, fallback_data_file):
    try:
        with open(param_file, 'r', encoding='latin1') as f: lines = f.readlines()
    except FileNotFoundError:
        with open(fallback_data_file, 'r', encoding='latin1') as f: lines = f.readlines()
    params = {}
    for line in lines:
        if line.startswith('#') or line.strip() == '': continue
        if '=' in line:
            key, val = map(str.strip, line.split('=', 1))
            try: params[key] = float(val)
            except ValueError: continue
    required = ['burgMag']
    for k in required:
        if k not in params: raise ValueError(f"Missing parameter: {k}")
    return params

def parse_dislocation_data(data_file):
    nodes, segments = {}, []
    with open(data_file, 'r', encoding='latin1') as f: lines = f.readlines()
    i, parsing = 0, False
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("nodalData"):
            parsing = True; i += 1; continue
        if parsing:
            if line.startswith('#') or line == "": i += 1; continue
            if re.match(r"\d+,\d+", line):
                parts = line.split()
                node_id, num_arms, constraint = int(parts[0].split(",")[1]), int(parts[4]), int(parts[5])
                x, y, z = map(float, parts[1:4])
                nodes[node_id] = {'coord': np.array([x, y, z]), 'constraint': constraint}
                for j in range(num_arms):
                    arm_line = lines[i + 1 + 2 * j].strip()
                    arm_node_id = int(arm_line.split()[0].split(",")[1])
                    burg = list(map(float, arm_line.split()[1:4]))
                    segments.append([node_id, arm_node_id, np.array(burg)])
                i += 1 + 2 * num_arms
            else: i += 1
        else: i += 1
    return nodes, segments

def export_vtp(filename, nodes, segments, point_data, cell_data):
    from pyevtk.hl import polyLinesToVTK
    point_list, points_per_line = [], []
    valid_segments = [s for s in segments if s[0] in nodes and s[1] in nodes]
    for n1, n2, _ in valid_segments:
        point_list.append(nodes[n1]['coord']); point_list.append(nodes[n2]['coord'])
        points_per_line.append(2)
    if not point_list: return
    points = np.array(point_list, dtype=np.float32)
    polyLinesToVTK(filename, x=np.ascontiguousarray(points[:, 0]), y=np.ascontiguousarray(points[:, 1]),
                   z=np.ascontiguousarray(points[:, 2]), pointsPerLine=np.array(points_per_line, dtype=np.int32),
                   pointData=point_data, cellData=cell_data)

def process_single_dislocation(i, args):
    base = f"rs{i:04d}"
    param_file = os.path.join(RESTART_DIR, base)
    data_file = os.path.join(RESTART_DIR, f"{base}.data")
    if not os.path.exists(data_file): return None
    out_path = os.path.join(PARADIS_VTK_DIR, f"{base}_dislocations")
    try:
        params = parse_params(param_file, data_file)
        burgMag = params.get('burgMag', 1.0)
        nodes, segments = parse_dislocation_data(data_file)
    except Exception as e:
        print(f"[ERROR] Could not parse {data_file}: {e}"); return None
    # for node_id in nodes: nodes[node_id]['coord'] *= burgMag
    valid_segments = [s for s in segments if s[0] in nodes and s[1] in nodes]
    if not valid_segments: return None
    burg_map, burg_id_counter = {}, 0
    burg_ids, burg_vectors_raw, constraints = [], [], []
    burg_tolerance_digits = 4
    for n1, n2, burg in valid_segments:
        burg_rounded = np.round(burg, burg_tolerance_digits)
        canonical_burg = burg_rounded.copy()
        for j in range(len(canonical_burg)):
            if canonical_burg[j] != 0:
                if canonical_burg[j] < 0: canonical_burg *= -1
                break
        burg_tuple = tuple(canonical_burg)
        if burg_tuple not in burg_map: burg_map[burg_tuple] = burg_id_counter; burg_id_counter += 1
        burg_ids.append(burg_map[burg_tuple])
        burg_vectors_raw.append(burg * burgMag)
        constraints.append(nodes[n1]['constraint'])
    burg_vectors_array = np.array(burg_vectors_raw, dtype=np.float32)
    cell_data = {
        "BurgersVectorID": np.array(burg_ids, dtype=np.int32), "Constraint": np.array(constraints, dtype=np.int32),
        "BurgersVector": (np.ascontiguousarray(burg_vectors_array[:, 0]), np.ascontiguousarray(burg_vectors_array[:, 1]),
                          np.ascontiguousarray(burg_vectors_array[:, 2]))}
    export_vtp(out_path, nodes, segments, {}, cell_data)
    return os.path.basename(out_path)

def write_paradis_pvd(pvd_path, basenames):
    with open(pvd_path, 'w') as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1">\n  <Collection>\n')
        for i, name in enumerate(sorted(basenames)):
            f.write(f'    <DataSet timestep="{i}" file="{name}.vtp"/>\n')
        f.write('  </Collection>\n</VTKFile>\n')
    print(f"Generated ParaDiS PVD file: {os.path.abspath(pvd_path)}")

# ==============================================================================
# ==                            Main Execution Block                        ==
# ==============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified script for Abaqus & ParaDiS VTK conversion.")
    subparsers = parser.add_subparsers(dest='mode', required=True, help='Select operation mode')

    full_run_parser = subparsers.add_parser('full_run', help='Run the full pipeline: Abaqus -> ParaDiS -> Zip.')
    full_run_parser.add_argument("--odbFile", type=str, required=True, help="Name of the Abaqus ODB file.")
    full_run_parser.add_argument("--no-parallel", action='store_true', help="Disable parallel processing.")
    full_run_parser.add_argument("--n-frames", type=int, default=None, help="Number of frames to sample from Abaqus and ParaDiS.")

    paradis_parser = subparsers.add_parser('paradis_only', help='Run only the ParaDiS to VTK conversion.')
    paradis_parser.add_argument("--paradis-frames", type=int, required=True, help="Number of frames to sample for ParaDiS.")
    paradis_parser.add_argument("--no-parallel", action='store_true', help="Disable parallel processing.")
    paradis_parser.add_argument("--n-frames", type=int, default=None, help="Override number of frames for ParaDiS only.")
    
    args = parser.parse_args()

    os.makedirs(PARADIS_VTK_DIR, exist_ok=True)
    os.makedirs(ZIP_OUTPUT_DIR, exist_ok=True)
    
    abaqus_output_dir = None
    paradis_output_dir = None

    worker_script_path = create_abaqus_worker_script()

    try:
        if args.mode == 'full_run':
            frame_count, abaqus_output_dir = abaqus_main(args, worker_script_path)
            if frame_count is None:
                print("[FATAL] Abaqus processing failed. Aborting."); sys.exit(1)
            
            paradis_output_dir = paradis_main(args, frame_count)
            if paradis_output_dir is None:
                print("[FATAL] ParaDiS processing failed. Aborting."); sys.exit(1)
                
        elif args.mode == 'paradis_only':
            paradis_output_dir = paradis_main(args)
            if paradis_output_dir is None:
                 print("[FATAL] ParaDiS processing failed. Aborting."); sys.exit(1)
        
        source_dirs_to_zip = [d for d in [abaqus_output_dir, paradis_output_dir] if d is not None]
        if source_dirs_to_zip:
            zip_filename = f"{os.path.splitext(args.odbFile)[0] if hasattr(args, 'odbFile') else 'paradis'}_vtk_results.zip"
            create_zip_archive(source_dirs_to_zip, zip_filename)

    finally:
        if os.path.exists(worker_script_path):
            os.remove(worker_script_path)
            print(f"Removed temporary worker script: {ABAQUS_WORKER_SCRIPT_NAME}")

    print("\n[SUCCESS] All tasks completed.")