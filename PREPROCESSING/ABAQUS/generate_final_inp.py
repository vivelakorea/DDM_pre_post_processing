import re
from collections import defaultdict
import os
import argparse

def generate_final_inp(input_file, final_output, ela_test=False):
    PART_NAME = 'Microbeam'
    Elset_col_num = 16

    # 1. Center alignment
    if not input_file.endswith('.inp'):
        input_file += '.inp'
    else:
        print("Input file doesn't need to end with '.inp'.")

    if not final_output.endswith('.inp'):
        final_output += '.inp'
    else:
        print("Output file doesn't need to end with '.inp'.")

    input_file = os.path.join("input", input_file)
    with open(input_file, "r") as f:
        lines = f.readlines()

    nodes = []
    in_node_block = False
    for line in lines:
        if line.lower().startswith("*node"):
            in_node_block = True
            continue
        if in_node_block:
            if line.startswith("*"):
                break
            parts = line.strip().split(",")
            if len(parts) >= 4:
                nid, x, y, z = int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
                nodes.append([nid, x, y, z])

    xs, ys, zs = zip(*[(n[1], n[2], n[3]) for n in nodes])
    cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)

    shifted_lines = []
    in_node_block = False
    for line in lines:
        if line.lower().startswith("*node"):
            shifted_lines.append(line)
            in_node_block = True
            continue
        if in_node_block:
            if line.startswith("*"):
                in_node_block = False
                shifted_lines.append(line)
                continue
            parts = line.strip().split(",")
            if len(parts) >= 4:
                nid, x, y, z = int(parts[0]), float(parts[1])-cx, float(parts[2])-cy, float(parts[3])-cz
                shifted_lines.append(f"{nid}, {x:.6e}, {y:.6e}, {z:.6e}\n")
            else:
                shifted_lines.append(line)
        else:
            shifted_lines.append(line)

    # 2. Extract surface elements (face-based)
    elements = {}
    current_section = None
    for line in shifted_lines:
        if line.lower().startswith("*element"):
            current_section = "element"
            continue
        if line.startswith("*"):
            current_section = None
        elif current_section == "element":
            parts = list(map(int, line.strip().split(",")))
            elements[parts[0]] = parts[1:]

    def get_faces(conn):
        return [tuple(sorted([conn[i] for i in face])) for face in [
            [0,1,2,3], [4,5,6,7], [0,1,5,4], [2,3,7,6], [0,3,7,4], [1,2,6,5]
        ]]

    face_count = defaultdict(list)
    for eid, conn in elements.items():
        for face in get_faces(conn):
            face_count[face].append(eid)

    surface_elements = sorted({eids[0] for face, eids in face_count.items() if len(eids) == 1})

    surfelem_lines = []
    surfelem_lines.append(f"*Elset, elset=Set-SurfaceEl, instance={PART_NAME}-1\n")
    for i in range(0, len(surface_elements), Elset_col_num):
        surfelem_lines.append("  " + ", ".join(str(eid) for eid in surface_elements[i:i+Elset_col_num]) + "\n")

    # 3. Define nodesets
    x_min = min(n[1] for n in nodes)
    x_max = max(n[1] for n in nodes)
    y_min = min(n[2] for n in nodes)
    y_max = max(n[2] for n in nodes)
    z_min = min(n[3] for n in nodes)
    z_max = max(n[3] for n in nodes)

    def find_node(x_cond=None, y_cond=None, z_cond=None):
        result = []
        for n in nodes:
            xok = True if x_cond is None else abs(n[1] - x_cond) < 1e-6
            yok = True if y_cond is None else abs(n[2] - y_cond) < 1e-6
            zok = True if z_cond is None else abs(n[3] - z_cond) < 1e-6
            if xok and yok and zok:
                result.append(n[0])
        return result

    coord_map = {
        'min': {'x': x_min, 'y': y_min, 'z': z_min},
        'max': {'x': x_max, 'y': y_max, 'z': z_max},
        None: None
    }

    choices = [None, 'min', 'max']
    nodeset_lines = []

    for x_choice in choices:
        for y_choice in choices:
            for z_choice in choices:
                if x_choice is None and y_choice is None and z_choice is None:
                    continue
                x_val = coord_map[x_choice]['x'] if x_choice else None
                y_val = coord_map[y_choice]['y'] if y_choice else None
                z_val = coord_map[z_choice]['z'] if z_choice else None
                matched_nodes = find_node(x_val, y_val, z_val)
                if not matched_nodes:
                    continue
                name_parts = []
                if x_choice: name_parts.append(f"X0" if x_choice == 'min' else "XL")
                if y_choice: name_parts.append(f"Y0" if y_choice == 'min' else "YL")
                if z_choice: name_parts.append(f"Z0" if z_choice == 'min' else "ZL")
                nset_name = "".join(name_parts)
                nodeset_lines.append(f"*Nset, nset={nset_name}, instance={PART_NAME}-1\n")
                for i in range(0, len(matched_nodes), Elset_col_num):
                    line = ", ".join(map(str, matched_nodes[i:i+Elset_col_num]))
                    nodeset_lines.append("  " + line + "\n")

    # 4. Define elset_all
    max_eid = max(elements.keys())
    all_elset_block = [
        "*Elset, elset=elset_all, generate\n",
        f"    1, {max_eid},     1\n"
    ]
    all_section_block = [
        "** Section: Section-1-elset_all\n",
        "*Solid Section, elset=elset_all, controls=EC-1, material=m1\n",
        ",\n"
    ]

    # 5. Replace text and merge
    shifted_lines = [line.replace("Part-1", PART_NAME) for line in shifted_lines]
    part_end_idx = next(i for i, line in enumerate(shifted_lines) if "*End Part" in line)
    pre_part = shifted_lines[:part_end_idx]
    post_part = shifted_lines[part_end_idx:]

    if ela_test:
        footer_files = ["footer1_ELA.inp", 
                        "footer2_ELA.inp", 
                        "footer3_BC0.inp", 
                        "footer4_ELA.inp", 
                        "footer5_BC1_ELA.inp", 
                        "footer6.inp"]
    else:
        footer_files = ["footer1_VUMAT.inp", 
                        "footer2_VUMAT.inp", 
                        "footer3_BC0.inp", 
                        "footer4_VUMAT.inp", 
                        "footer5_BC1_VUMAT.inp", 
                        "footer6.inp"]

    footers = []
    for footer_file in footer_files:
        path = os.path.join(".","input", "footer", footer_file)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        
        with open(path, "r") as f:
            footers.append(f.readlines())


    final = (
        pre_part +
        all_elset_block +
        all_section_block +
        post_part +
        surfelem_lines +
        nodeset_lines
    )
    for footer in footers:
        final += footer
        final += ["\n"]

    final_output = os.path.join("output", "Job-" + final_output)
    if not os.path.exists(os.path.dirname(final_output)):
        os.makedirs(os.path.dirname(final_output))

    with open(final_output, "w") as f:
        f.writelines(final)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final inp file.")
    parser.add_argument("--input", required=True, help="Input .inp file")
    parser.add_argument("--output", required=True, help="Output .inp file")
    parser.add_argument("--ela_test", action='store_true', help="Generate for ELA test case", default=False)
    args = parser.parse_args()
    generate_final_inp(args.input, args.output, ela_test=args.ela_test)
