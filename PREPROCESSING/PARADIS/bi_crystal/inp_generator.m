%% ParaDiS Bicrystal Frank-Read Source Network Generator
%
% Generates ParaDiS input files (.data, .ctrl) where all geometric data in
% the .data file is in units of Burgers vector magnitude (b).
% Creates a bicrystal structure with a grain boundary on the XZ plane (y=0).

clc;
clear;
close all;
rng(42); % for reproducibility

%% ========================================================================
%  1. USER INPUT PARAMETERS
%  ========================================================================


config.test_paradis_only = false;


% --- Script Execution Control ---
config.num_files_to_generate = 1;
config.enable_plotting = true;


% --- Crystal Structure & Material ---
config.crystalStructure = 'HCP';
config.mobilityLaw = 'HCP_linear';

% --- Bicrystal Rotation (Bunge Euler Angles) ---
bunge_phi1 = 45.0; bunge_PHI = 54.74; bunge_phi2 = 0.0; % In DEGREES
bunge_phi1 = 0.0; bunge_PHI = 90.0; bunge_phi2 = 0.0; % For HCP soft grain
config.R = bunge_euler_to_matrix(deg2rad(bunge_phi1), deg2rad(bunge_PHI), deg2rad(bunge_phi2));
config.R = (config.R)';

% --- Geometry and Density ---
% Box size is in units of b.
config.boxSize_b = [3500,3500,7000];
% NOTE: 1e5 is very low. Typical values are 1e12 to 1e14.
config.targetDensity = 1e13; % Target dislocation density in [m^-2]

% --- Frank-Read Source Parameters ---
config.frs_length_alpha = 1.00; % mean_frs_length_m = config.frs_length_alpha / sqrt(config.targetDensity);
config.frs_length_sigma = 0.00; % frs_length = mean_frs_length_b + std_dev_b * randn();

config.maxSeg_b = 391.236307;
config.minSeg_b = 78.247261;

%
config.name = 'D1750_Ti';


% --- Material and Mobility Parameters ---
config.burgMag = 2.556e-10; % Meters
config.shearModulus = 4.8e10; % Pa
config.pois = 0.34;

if strcmp(config.crystalStructure, "HCP")
    % Be

    config.rc=3.00000; 

    % HCP specific parameters
    % cOVERa, HCPEcoreA, HCPEcoreC, HCPEcoreCpA, 
    
    config.cOVERa = 1.58; % c/a ratio
    config.HCPEcoreA = 18.06977e9; % Ecore energy for A
    config.HCPEcoreC = 21.70031e9; % Ecore energy for C
    config.HCPEcoreCpA = 7.18520e9; % Ecore energy for CpA
    
    % HCP_A_Basal_EdgeDrag, HCP_A_Basal_ScrewDrag,
    % HCP_A_Prismatic_EdgeDrag, HCP_A_Prismatic_ScrewDrag,
    % HCP_A_1stPyramidal_EdgeDrag, HCP_A_1stPyramidal_ScrewDrag,
    % HCP_A_2ndPyramidal_EdgeDrag, HCP_A_2ndPyramidal_ScrewDrag,

    config.HCP_A_Basal_EdgeDrag = 1.18434387898773e-04; % Basal edge drag
    config.HCP_A_Basal_ScrewDrag = 1.65421659553579e-04; % Basal screw drag
    config.HCP_A_Prismatic_EdgeDrag = 5.69513923115224e-05; % Prismatic edge drag
    config.HCP_A_Prismatic_ScrewDrag = 7.15096201802962e-05; % Prismatic screw drag
    config.HCP_A_1stPyramidal_EdgeDrag = 7.28268834977848e-05; % 1st Pyramidal edge drag
    config.HCP_A_1stPyramidal_ScrewDrag = 1.05154866364405e-04; % 1st Pyramidal screw drag
    config.HCP_A_2ndPyramidal_EdgeDrag = 7.28268834977848e-05; % 2nd Pyramidal edge drag
    config.HCP_A_2ndPyramidal_ScrewDrag = 1.05154866364405e-04; % 2nd Pyramidal screw drag

    % From polyxtal literature, values are assumed 20 times basal edge
    % HCP_CpA_Prismatic_EdgeDrag, HCP_CpA_Prismatic_ScrewDrag,
    % HCP_CpA_1stPyramidal_EdgeDrag, HCP_CpA_1stPyramidal_ScrewDrag,
    % HCP_CpA_2ndPyramidal_EdgeDrag, HCP_CpA_2ndPyramidal_ScrewDrag,
    % HCP_C_Prismatic_EdgeDrag, HCP_C_Prismatic_ScrewDrag,

    
    config.HCP_CpA_Prismatic_EdgeDrag = config.HCP_A_Basal_EdgeDrag * 20; % CpA Prismatic edge drag
    config.HCP_CpA_Prismatic_ScrewDrag = config.HCP_A_Basal_EdgeDrag * 20; % CpA Prismatic screw drag
    config.HCP_CpA_1stPyramidal_EdgeDrag = config.HCP_A_Basal_EdgeDrag * 20; % CpA 1st Pyramidal edge drag
    config.HCP_CpA_1stPyramidal_ScrewDrag = config.HCP_A_Basal_EdgeDrag * 20; % CpA 1st Pyramidal screw drag
    config.HCP_CpA_2ndPyramidal_EdgeDrag = config.HCP_A_Basal_EdgeDrag * 20; % CpA 2nd Pyramidal edge drag
    config.HCP_CpA_2ndPyramidal_ScrewDrag = config.HCP_A_Basal_EdgeDrag * 20; % CpA 2nd Pyramidal screw drag
    config.HCP_C_Prismatic_EdgeDrag = config.HCP_A_Basal_EdgeDrag * 20; % C Prismatic edge drag
    config.HCP_C_Prismatic_ScrewDrag = config.HCP_A_Basal_EdgeDrag * 20; % C Prismatic screw drag

    % 10000 * pyramidal
    % HCP_Sessile_EdgeDrag, HCP_Sessile_ScrewDrag,
    
    config.HCP_Sessile_EdgeDrag = 10000 * config.HCP_A_1stPyramidal_EdgeDrag; % Sessile edge drag
    config.HCP_Sessile_ScrewDrag = 10000 * config.HCP_A_1stPyramidal_ScrewDrag; % Sessile screw drag

    % ~10 times lowest
    % HCP_LineDrag

    config.HCP_LineDrag = 10 * config.HCP_A_Basal_EdgeDrag; % Line drag

else
    config.rc = 0.100000; 

    % FCC/BCC specific parameters
    config.MobScrew = 1.0e+05; % Mobility for screw dislocations
    config.MobEdge = 1.0e+05; % Mobility for edge dislocations
    config.MobGB = 5.0; % Grain boundary mobility
end
config.MobScrew=1.0e+05; 
config.MobEdge=1.0e+05; 
config.MobGB=5.0;


% ... (rest of the config is the same) ...


config.numXdoms=1; config.numYdoms=1; config.numZdoms=1; config.numXcells=4;
config.numYcells=4; config.numZcells=4; config.xBoundType=0; config.yBoundType=0;
config.zBoundType=0; config.DLBfreq=0; 
config.timestepIntegrator='explicit'; config.maxDT=1e-12; config.rTol=15.6495;
config.remeshRule=2; config.splitMultiNodeFreq=1;
config.fmEnabled=0; config.fmMPOrder=2; config.fmTaylorOrder=5;
config.fmCorrectionTbl="inputs/fm-ctab.Ta.600K.0GPa.m2.t5.data";
config.Rijmfile="inputs/Rijm.cube.out"; config.RijmPBCfile="inputs/RijmPBC.cube.out";
config.winDefaultsFile="inputs/paradis.xdefaults.thinfilm"; config.useLabFrame=1;
config.TempK=300; config.pressure=0.0; config.loadType=1;
config.appliedStress=[0;0;0;0;0;0]; config.indxErate=0; config.vacancyConc=-1.0;
config.vacancyConcEquilibrium=-1.0; config.MobClimb=1e-08;
config.MobGlide=1e-08; config.includeInertia=0; config.massDensity=-1.0;
config.armfile=1; config.armfilefreq=1; config.savecn=1; config.savecnfreq=1;
config.savecncounter=0; config.saveprop=1; config.savepropfreq=1;
config.atomeye=1; config.atomeyefreq=100; config.atomeyesegradius=100;
config.writeFlux=1; config.writeFluxFullDecomp=1;
config.writeFluxFullDecompTotals=1; config.writeFluxSimpleTotals=1;
config.writeFluxFreq=1; config.FEM_DD_ForceSegTol=2000;
config.enforceGlidePlanes=1; config.enableCrossSlip=0;
config.TensionFactor=1.0; config.elasticinteraction=1;


if config.test_paradis_only
    config.maxstep=100;
    config.edotdir=[1.0;0.0;0.0];
    config.eRate=1e+10;
else
    config.maxstep=1;
    config.edotdir=[1.0;0.0;0.0];
    config.eRate=0;
end


config.output_folder = 'output';

%% ========================================================================
%  2. SCRIPT EXECUTION LOGIC
%  ========================================================================

% --- Derived Parameters ---
box_b = config.boxSize_b;
targetDensity_b2 = config.targetDensity * config.burgMag^2;
volume_b3 = prod(box_b);

% --- Get Slip Systems ---
slipSystems = get_slip_systems(config);

% --- Main loop ---
for i_file = 1:config.num_files_to_generate
    fprintf('--- Generating bicrystal file set %d of %d ---\n', i_file, config.num_files_to_generate);
    rn = []; links = []; total_length_b = 0; current_density_b2 = 0;
    while current_density_b2 < targetDensity_b2
        center_pos = (rand(1, 3) - 0.5) .* box_b;
        rand_idx = randi(length(slipSystems));
        b_base = slipSystems(rand_idx).b; n_base = slipSystems(rand_idx).n;
        if center_pos(2) < 0, b = b_base; n = n_base;
        else, b = (config.R * b_base')'; n = (config.R * n_base')'; end
        if rand() <= 0.5, lineDirection = b; else, lineDirection = cross(b, n); end
        lineDirection = lineDirection / norm(lineDirection);
        
        mean_frs_length_m = config.frs_length_alpha / sqrt(config.targetDensity);
        mean_frs_length_b = mean_frs_length_m / config.burgMag;
        std_dev_b = config.frs_length_sigma * mean_frs_length_b;
        frs_length = mean_frs_length_b + std_dev_b * randn();

        % +++ BUG FIX: ADDED SAFETY CHECK FOR FRS LENGTH +++
        % Cap the FRS length to prevent it from exceeding the box size.
        max_possible_length = min(box_b) * 0.8; % Max length is 80% of the smallest box dimension
        if frs_length > max_possible_length
            frs_length = max_possible_length * (0.8 + 0.4*rand()); % Use a random fraction of the max length
        end
        % +++ END OF BUG FIX +++

        if frs_length <= 0, continue; end
        
        p1 = center_pos - lineDirection * frs_length / 2;
        p3 = center_pos + lineDirection * frs_length / 2;
        if any(abs(p1) > box_b / 2) || any(abs(p3) > box_b / 2), continue; end
        if (p1(2) * p3(2)) < 0, continue; end
        
        p2 = center_pos;
        node_id_start = size(rn, 1);
        rn = [rn; [p1, 7]; [p2, 0]; [p3, 7]];
        links = [links; [node_id_start+1, node_id_start+2, b, n]; [node_id_start+2, node_id_start+3, b, n]];
        total_length_b = total_length_b + frs_length;
        current_density_b2 = total_length_b / volume_b3;
    end
    fprintf('Target density reached. Final Density: %4.4e [/m^2]\n', current_density_b2 / (config.burgMag^2));
    % base_filename = sprintf('Bicrystal_%s_run_%d', config.crystalStructure, i_file);
    base_filename = sprintf('%s_%d', config.name, i_file);
    
    % 출력 폴더 생성
    if ~exist(config.output_folder, 'dir')
        mkdir(config.output_folder);
    end
    
    % 결과 폴더 경로 생성 및 생성
    results_folder = fullfile(config.output_folder, sprintf('%s_results', base_filename));
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
    
    % 파일 경로 지정
    ctrl_path = fullfile(results_folder, [base_filename, '.ctrl']);
    data_path = fullfile(results_folder, [base_filename, '.data']);
    
    % 디버그 출력
    disp(ctrl_path);
    
    % 파일 생성 함수 호출
    write_paradis_control_file(ctrl_path, config, base_filename);
    write_paradis_data_file(config, data_path, rn, links, box_b);
    end
    
    if config.enable_plotting, plot_bicrystal(rn, links, box_b, config); end


%% ========================================================================
%  3. LOCAL FUNCTIONS
%  ========================================================================

function R = bunge_euler_to_matrix(phi1, PHI, phi2)
    c1=cos(phi1); s1=sin(phi1); c_=cos(PHI); s_=sin(PHI); c2=cos(phi2); s2=sin(phi2);
    R = [c1*c2-s1*s2*c_, -c1*s2-s1*c2*c_, s1*s_; s1*c2+c1*s2*c_, -s1*s2+c1*c2*c_, -c1*s_; s2*s_, c2*s_, c_];
end

function systems = get_slip_systems(config)
    crystal_structure = config.crystalStructure;
    
    % Initialize the struct array
    systems = struct('b', {}, 'n', {});

    if strcmpi(crystal_structure, 'HCP')
        c_a = config.cOVERa;
        a = 1.0; 
        c = c_a * a;

        % --- Base Cartesian Vectors for burgers and normals ---
        % a-type vectors <11-20>/3 type
        b1 = [a, 0, 0];
        b2 = [-a/2, a*sqrt(3)/2, 0];
        b3 = [-a/2, -a*sqrt(3)/2, 0];
        
        % <c+a>-type vectors <11-23>/3 type
        b4 = [a, 0, c];
        b5 = [-a/2,  a*sqrt(3)/2, c];
        b6 = [-a/2, -a*sqrt(3)/2, c];
        b7 = [a, 0, -c];
        b8 = [-a/2,  a*sqrt(3)/2, -c];
        b9 = [-a/2, -a*sqrt(3)/2, -c];

        % <c>-type Burgers vector [0001]
        b10 = [0, 0, c];

        % --- Slip Plane Normals ---
        % Basal (0001)
        n_basal = [0, 0, 1];
        % Prismatic {10-10}
        n_pr1 = [0, 1, 0]; n_pr2 = [sqrt(3)/2, -1/2, 0]; n_pr3 = [-sqrt(3)/2, -1/2, 0];
        % 1st Pyramidal {10-11}
        n_pyI1 = [0, -2*c/a, sqrt(3)]; n_pyI2 = [0, 2*c/a, sqrt(3)];
        n_pyI3 = [-sqrt(3)*c/a, -c/a, -sqrt(3)]; n_pyI4 = [sqrt(3)*c/a, c/a, -sqrt(3)];
        n_pyI5 = [sqrt(3)*c/a, -c/a, sqrt(3)]; n_pyI6 = [-sqrt(3)*c/a, c/a, sqrt(3)];
        % 2nd Pyramidal {11-22}
        n_pII1 = [a, 0, -c]; n_pII2 = [-a/2, a*sqrt(3)/2, -c]; n_pII3 = [-a/2, -a*sqrt(3)/2, -c];
        n_pII4 = [-a, 0, c]; n_pII5 = [a/2, -a*sqrt(3)/2, c]; n_pII6 = [a/2, a*sqrt(3)/2, c];

        % =====================================================================
        % Add systems with orthogonality check
        % =====================================================================
        
        % Table 1: Burgers vectors of type <a> (Total: 12 systems)
        % Basal (3 systems)
        add_system_if_orthogonal(b1, n_basal);
        add_system_if_orthogonal(b2, n_basal);
        add_system_if_orthogonal(b3, n_basal);

        % Prismatic (3 systems)
        add_system_if_orthogonal(b1, n_pr1);
        add_system_if_orthogonal(b3, n_pr2);
        add_system_if_orthogonal(b2, n_pr3);

        % 1st Pyramidal (6 systems)
        add_system_if_orthogonal(b1, n_pyI1); add_system_if_orthogonal(b1, n_pyI2);
        add_system_if_orthogonal(b2, n_pyI3); add_system_if_orthogonal(b2, n_pyI4);
        add_system_if_orthogonal(b3, n_pyI5); add_system_if_orthogonal(b3, n_pyI6);

        % Table 2: Burgers vectors of type <c+a> (Total: 24 systems)
        % Prismatic (6 systems)
        add_system_if_orthogonal(b4, n_pr1); add_system_if_orthogonal(b7, n_pr1);
        add_system_if_orthogonal(b5, n_pr3); add_system_if_orthogonal(b8, n_pr3);
        add_system_if_orthogonal(b6, n_pr2); add_system_if_orthogonal(b9, n_pr2);
        
        % 1st Pyramidal (12 systems)
        add_system_if_orthogonal(b5, n_pyI1); add_system_if_orthogonal(b6, n_pyI1);
        add_system_if_orthogonal(b8, n_pyI2); add_system_if_orthogonal(b9, n_pyI2);
        add_system_if_orthogonal(b6, n_pyI3); add_system_if_orthogonal(b7, n_pyI3);
        add_system_if_orthogonal(b5, n_pyI4); add_system_if_orthogonal(b4, n_pyI4);
        add_system_if_orthogonal(b4, n_pyI5); add_system_if_orthogonal(b8, n_pyI5);
        add_system_if_orthogonal(b7, n_pyI6); add_system_if_orthogonal(b9, n_pyI6);

        % 2nd Pyramidal (6 systems)
        add_system_if_orthogonal(b5, n_pII1);
        add_system_if_orthogonal(b6, n_pII2);
        add_system_if_orthogonal(b4, n_pII3);
        add_system_if_orthogonal(b8, n_pII4);
        add_system_if_orthogonal(b9, n_pII5);
        add_system_if_orthogonal(b7, n_pII6);
        
        % Table 3: Burgers vector of type <c> (Total: 3 systems)
        add_system_if_orthogonal(b10, n_pr1);
        add_system_if_orthogonal(b10, n_pr2);
        add_system_if_orthogonal(b10, n_pr3);

    else % FCC or BCC Logic (Unchanged)
        if strcmpi(crystal_structure, 'FCC')
            planes = [[1 1 1]; [1 1 -1]; [1 -1 1]; [-1 1 1]];
            dirs = [[1 -1 0]; [1 1 0]; [1 0 -1]; [1 0 1]; [0 1 -1]; [0 1 1]];
        elseif strcmpi(crystal_structure, 'BCC')
            planes = [[1 1 0]; [1 -1 0]; [1 0 1]; [1 0 -1]; [0 1 1]; [0 1 -1]];
            dirs = [[1 1 1]; [1 1 -1]; [1 -1 1]; [-1 1 1]];
        else
             error('Unsupported crystal structure: %s', crystal_structure);
        end
        for p = 1:size(planes,1)
            for d = 1:size(dirs,1)
                if abs(dot(planes(p,:), dirs(d,:))) < 1e-9
                    systems(end+1).b = dirs(d,:) / norm(dirs(d,:));
                    systems(end).n = planes(p,:) / norm(planes(p,:));
                end
            end
        end
    end
    
    % --- Nested helper function for orthogonality check and addition ---
    function add_system_if_orthogonal(b_vec, n_vec)
        if abs(dot(b_vec, n_vec)) > 1e-9
            fprintf('WARNING: Non-orthogonal pair skipped. dot(b,n) = %.4f\n', dot(b_vec, n_vec));
            fprintf('  b = [%.4f, %.4f, %.4f]\n', b_vec);
            fprintf('  n = [%.4f, %.4f, %.4f]\n', n_vec);
        else
            systems(end+1).b = b_vec;
            systems(end).n = n_vec / norm(n_vec);
        end
    end
end
function write_paradis_data_file(config, fname, rn, links, box_b)

    if strcmp(config.crystalStructure, 'HCP')
        Nnodes = size(rn, 1);
        LINKMAX = size(links, 1);
        half_box_b = box_b / 2;
        
        list = cell(Nnodes, 1);
        for j = 1:LINKMAX
            n0 = links(j, 1);
            n1 = links(j, 2);
            list{n0} = [list{n0}, j];
            list{n1} = [list{n1}, j];
        end
        
        fid = fopen(fname, 'w');
        if fid == -1, error('Cannot open file: %s', fname); end

        if config.test_paradis_only
            max_half_box_size = max(half_box_b)*4;
        else
            max_half_box_size = max(half_box_b);
        end
    
        
        % Corrected Header for v5 format
        fprintf(fid, 'dataFileVersion =   5\n');
        fprintf(fid, 'numFileSegments =   1\n');
        fprintf(fid, 'minCoordinates = [\n  %e\n  %e\n  %e\n  ]\n', -max_half_box_size, -max_half_box_size, -max_half_box_size);
        fprintf(fid, 'maxCoordinates = [\n  %e\n  %e\n  %e\n  ]\n',  max_half_box_size, max_half_box_size, max_half_box_size);
        fprintf(fid, 'nodeCount =   %d\n', Nnodes);
        fprintf(fid, 'dataDecompType =   2\n');
        fprintf(fid, 'dataDecompGeometry = [\n  1\n  1\n  1\n  ]\n\n');
        
        fprintf(fid, '#\n#  END OF DATA FILE PARAMETERS\n#\n\n');
        
        % Corrected Domain Decomposition format
        fprintf(fid, 'domainDecomposition = \n');
        fprintf(fid, '# Dom_ID  Minimum XYZ bounds   Maximum XYZ bounds\n');
        fprintf(fid, '  0   %.4f  %.4f  %.4f    %.4f   %.4f   %.4f\n', ...
            -max_half_box_size, -max_half_box_size, -max_half_box_size, ...
             max_half_box_size,  max_half_box_size,  max_half_box_size);
             
        fprintf(fid, 'nodalData = \n');
        fprintf(fid, '#  Primary lines: node_tag, x, y, z, num_arms, constraint\n');
        fprintf(fid, '#  Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n');
    
        for i = 1:Nnodes
            numNbrs = length(list{i});
            constraint = rn(i, 4);
            node_pos_b = rn(i, 1:3);
            
            % Primary Line
            fprintf(fid, ' 9,%d %.4f %.4f %.4f %d %d\n', ...
                     i-1, node_pos_b(1), node_pos_b(2), node_pos_b(3), numNbrs, constraint);
            
            for j = 1:numNbrs
                link_idx = list{i}(j);
                n0 = links(link_idx, 1);
                n1 = links(link_idx, 2);
                bv = links(link_idx, 3:5);
                nv = links(link_idx, 6:8);
                
                if n0 == i
                    neighbor_node_idx = n1;
                else
                    neighbor_node_idx = n0;
                    bv = -bv; % Reverse Burgers vector for the opposite direction
                end
                
                % Secondary Lines (matching Be.txt format)
                fprintf(fid, '   9,%d %e %e %e\n', ...
                         neighbor_node_idx-1, bv(1), bv(2), bv(3));
                fprintf(fid, '       %e %e %e\n', ...
                         nv(1), nv(2), nv(3));
            end
        end
        
        fclose(fid);
        fprintf('ParaDiS data file (v5 format) written to: %s\n', fname);
    else

        Nnodes = size(rn, 1); LINKMAX = size(links, 1);
        half_box_b = box_b / 2;
        list = cell(Nnodes, 1);
        for j=1:LINKMAX, n0=links(j,1); n1=links(j,2); list{n0}=[list{n0},j]; list{n1}=[list{n1},j]; end
        fid = fopen(fname, 'w'); if fid == -1, error('Cannot open file: %s', fname); end
        
        fprintf(fid, '#\n#\tParaDiS nodal data file (by write_loop_data.m) \n#\n \ndataFileVersion = 4\n');
        fprintf(fid, 'numFileSegments = 1\n');
    
        if config.test_paradis_only
            max_half_box_size = max(half_box_b)*10;
        else
            max_half_box_size = max(half_box_b);
        end
    
        fprintf(fid, 'minCoordinates = [\n %e\n %e\n %e\n ]\n', -max_half_box_size, -max_half_box_size, -max_half_box_size);
        fprintf(fid, 'maxCoordinates = [\n %e\n %e\n %e\n ]\n',  max_half_box_size, max_half_box_size, max_half_box_size);
        fprintf(fid, 'nodeCount = %d\n', Nnodes);
        fprintf(fid, 'dataDecompType = 1\n');
        fprintf(fid, 'dataDecompGeometry = [\n 1\n 1\n 1\n ]\n');
        fprintf(fid, 'domainDecomposition = \n %e\n     %e\n         %e\n         %e\n     %e\n %e\n \n',...
            -max_half_box_size, -max_half_box_size, -max_half_box_size, ...
            max_half_box_size, max_half_box_size, max_half_box_size);
        fprintf(fid, '#\n#\tPrimary lines: node_tag, x, y, z, num_arms, constraint\n');
        fprintf(fid, '#\tSecondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n#\n');
        fprintf(fid, '#       length in unit of burgMag\n \nnodalData =\n');
    
        for i = 1:Nnodes
            numNbrs = length(list{i}); constraint = rn(i, 4);
            node_pos_b = rn(i, 1:3);
            fprintf(fid, '     %d,%-12d% 21.14e % 21.14e % 21.14e %d %d\n', ...
                     0, i-1, node_pos_b(1), node_pos_b(2), node_pos_b(3), numNbrs, constraint);
            for j = 1:numNbrs
                link_idx = list{i}(j); n0 = links(link_idx,1); n1 = links(link_idx,2);
                bv = links(link_idx, 3:5); nv = links(link_idx, 6:8);
                if n0 == i, neighbor_node_idx = n1; else, neighbor_node_idx = n0; bv = -bv; end
                fprintf(fid, '           %d,%-12d %-9.6f %-9.6f %-9.6f\n', ...
                         0, neighbor_node_idx-1, bv(1), bv(2), bv(3));
                fprintf(fid, '                 %-9.6f %-9.6f %-9.6f\n', ...
                         nv(1), nv(2), nv(3));
            end
        end
        fclose(fid); fprintf('ParaDiS data file written to: %s\n', fname);



    end
    


end

function write_paradis_control_file(fname, config, base_filename)
    fid = fopen(fname, 'w'); if fid == -1, error('Cannot open file: %s', fname); end
    
    fprintf(fid, '########################################\n###                                  ###\n');
    fprintf(fid, '###  ParaDiS control parameter file  ###\n###                                  ###\n');
    fprintf(fid, '########################################\n\n');
    fprintf(fid, 'dirname    = %s_results\n\n', ['tests/', base_filename]);


    fprintf(fid, '# Simulation cell and processor setup\n\n');
    fprintf(fid, 'numXdoms            = %d\n', config.numXdoms);
    fprintf(fid, 'numYdoms            = %d\n', config.numYdoms);
    fprintf(fid, 'numZdoms            = %d\n', config.numZdoms);
    fprintf(fid, 'numXcells           = %d\n', config.numXcells);
    fprintf(fid, 'numYcells           = %d\n', config.numYcells);
    fprintf(fid, 'numZcells           = %d\n', config.numZcells);
    fprintf(fid, 'xBoundType          = %d\n', config.xBoundType);
    fprintf(fid, 'yBoundType          = %d\n', config.yBoundType);
    fprintf(fid, 'zBoundType          = %d\n', config.zBoundType);
    fprintf(fid, 'DLBfreq             = %d\n\n', config.DLBfreq);


    fprintf(fid, '# Simulation time and timestepping controls\n\n');
    fprintf(fid, 'maxstep             = %d\n', config.maxstep);
    fprintf(fid, 'timestepIntegrator  = "%s"\n', config.timestepIntegrator);
    fprintf(fid, 'maxDT               = %e\n\n', config.maxDT);
    fprintf(fid, 'rTol                = %f\n', config.rTol);
    fprintf(fid, 'rc                  = %f\n\n', config.rc);


    fprintf(fid, '# Discretization and topological change controls\n\n');
    fprintf(fid, 'maxSeg              = %f\n', config.maxSeg_b);
    fprintf(fid, 'minSeg              = %f\n', config.minSeg_b);
    fprintf(fid, 'remeshRule          = %d\n', config.remeshRule);
    fprintf(fid, 'splitMultiNodeFreq  = %d\n\n', config.splitMultiNodeFreq);


    fprintf(fid, '# Fast multipole method specs.\n\n');
    fprintf(fid, 'fmEnabled           = %d\n', config.fmEnabled);
    fprintf(fid, 'fmMPOrder           = %d\n', config.fmMPOrder);
    fprintf(fid, 'fmTaylorOrder       = %d\n', config.fmTaylorOrder);
    fprintf(fid, 'fmCorrectionTbl     = "%s"\n\n', config.fmCorrectionTbl);


    fprintf(fid, '# Tables for non-FMM far-field force calcs.\n\n');
    fprintf(fid, 'Rijmfile            = "%s"\n', config.Rijmfile);
    fprintf(fid, 'RijmPBCfile         = "%s"\n', config.RijmPBCfile);
    fprintf(fid, 'winDefaultsFile     = "%s"\n\n', config.winDefaultsFile);


    fprintf(fid, '# Lab Frame\n\n');
    fprintf(fid, 'useLabFrame         = %d\n\n', config.useLabFrame);
    if config.useLabFrame == 1 && isfield(config, 'R')
        fprintf(fid, 'rotationMatrix = [\n');
        fprintf(fid, '  %-10f  %-10f  %-10f  \n', config.R(1,:));
        fprintf(fid, '  %-10f  %-10f  %-10f  \n', config.R(2,:));
        fprintf(fid, '  %-10f  %-10f  %-10f  \n', config.R(3,:));
        fprintf(fid, '  ]\n');
    end
    fprintf(fid, '\n');

    
    fprintf(fid, '# Loading conditions\n#\n');
    fprintf(fid, 'TempK               = %f\n', config.TempK);
    fprintf(fid, 'pressure            = %f\n', config.pressure);
    fprintf(fid, 'loadType            = %d \n\n', config.loadType);
    fprintf(fid, 'appliedStress = [\n  %f\n  %f\n  %f\n  %f\n  %f\n  %f\n  ]\n\n', config.appliedStress);
    fprintf(fid, 'eRate               = %e\n', config.eRate);
    fprintf(fid, 'indxErate           = %d\n', config.indxErate);
    fprintf(fid, 'edotdir = [ \n  %.1f\n  %.1f\n  %.1f\n  ]\n\n', config.edotdir);


    fprintf(fid, '# Material and mobility parameters\n\n');
    fprintf(fid, 'mobilityLaw     = "%s"\n', config.mobilityLaw);
    fprintf(fid, 'vacancyConc         = %f\n', config.vacancyConc);
    fprintf(fid, 'vacancyConcEquilibrium = %f\n\n', config.vacancyConcEquilibrium);
    fprintf(fid, 'shearModulus        = %e\n', config.shearModulus);
    fprintf(fid, 'pois                = %f\n', config.pois);
    fprintf(fid, 'burgMag             = %e\n\n', config.burgMag);

    if strcmp(config.crystalStructure, "HCP")
        % cOVERa, HCPEcoreA, HCPEcoreC, HCPEcoreCpA, 
        % HCP_A_Basal_EdgeDrag, HCP_A_Basal_ScrewDrag,
        % HCP_A_Prismatic_EdgeDrag, HCP_A_Prismatic_ScrewDrag,
        % HCP_A_1stPyramidal_EdgeDrag, HCP_A_1stPyramidal_ScrewDrag,
        % HCP_A_2ndPyramidal_EdgeDrag, HCP_A_2ndPyramidal_ScrewDrag,

        fprintf(fid, 'cOVERa              = %f\n', config.cOVERa);
        fprintf(fid, 'HCPEcoreA           = %e\n', config.HCPEcoreA);
        fprintf(fid, 'HCPEcoreC           = %e\n', config.HCPEcoreC);
        fprintf(fid, 'HCPEcoreCpA         = %e\n\n', config.HCPEcoreCpA);
        fprintf(fid, 'HCP_A_Basal_EdgeDrag = %e\n', config.HCP_A_Basal_EdgeDrag);
        fprintf(fid, 'HCP_A_Basal_ScrewDrag = %e\n', config.HCP_A_Basal_ScrewDrag);
        fprintf(fid, 'HCP_A_Prismatic_EdgeDrag = %e\n', config.HCP_A_Prismatic_EdgeDrag);
        fprintf(fid, 'HCP_A_Prismatic_ScrewDrag = %e\n', config.HCP_A_Prismatic_ScrewDrag);
        fprintf(fid, 'HCP_A_1stPyramidal_EdgeDrag = %e\n', config.HCP_A_1stPyramidal_EdgeDrag);
        fprintf(fid, 'HCP_A_1stPyramidal_ScrewDrag = %e\n', config.HCP_A_1stPyramidal_ScrewDrag);
        fprintf(fid, 'HCP_A_2ndPyramidal_EdgeDrag = %e\n', config.HCP_A_2ndPyramidal_EdgeDrag);
        fprintf(fid, 'HCP_A_2ndPyramidal_ScrewDrag = %e\n', config.HCP_A_2ndPyramidal_ScrewDrag);


        % From polyxtal literature, values are assumed 20 times basal edge
        % HCP_CpA_Prismatic_EdgeDrag, HCP_CpA_Prismatic_ScrewDrag,
        % HCP_CpA_1stPyramidal_EdgeDrag, HCP_CpA_1stPyramidal_ScrewDrag,
        % HCP_CpA_2ndPyramidal_EdgeDrag, HCP_CpA_2ndPyramidal_ScrewDrag,
        % HCP_C_Prismatic_EdgeDrag, HCP_C_Prismatic_ScrewDrag,

        fprintf(fid, 'HCP_CpA_Prismatic_EdgeDrag = %e\n', config.HCP_CpA_Prismatic_EdgeDrag);
        fprintf(fid, 'HCP_CpA_Prismatic_ScrewDrag = %e\n', config.HCP_CpA_Prismatic_ScrewDrag);
        fprintf(fid, 'HCP_CpA_1stPyramidal_EdgeDrag = %e\n', config.HCP_CpA_1stPyramidal_EdgeDrag);
        fprintf(fid, 'HCP_CpA_1stPyramidal_ScrewDrag = %e\n', config.HCP_CpA_1stPyramidal_ScrewDrag);
        fprintf(fid, 'HCP_CpA_2ndPyramidal_EdgeDrag = %e\n', config.HCP_CpA_2ndPyramidal_EdgeDrag);
        fprintf(fid, 'HCP_CpA_2ndPyramidal_ScrewDrag = %e\n', config.HCP_CpA_2ndPyramidal_ScrewDrag);
        fprintf(fid, 'HCP_C_Prismatic_EdgeDrag = %e\n', config.HCP_C_Prismatic_EdgeDrag);
        fprintf(fid, 'HCP_C_Prismatic_ScrewDrag = %e\n\n', config.HCP_C_Prismatic_ScrewDrag);


        % 10000 * pyramidal
        % HCP_Sessile_EdgeDrag, HCP_Sessile_ScrewDrag,

        fprintf(fid, 'HCP_Sessile_EdgeDrag = %e\n', config.HCP_Sessile_EdgeDrag);
        fprintf(fid, 'HCP_Sessile_ScrewDrag = %e\n\n', config.HCP_Sessile_ScrewDrag);

        % ~10 times lowest
        % HCP_LineDrag
        fprintf(fid, 'HCP_LineDrag = %e\n\n', config.HCP_LineDrag);


    else
        fprintf(fid, 'MobScrew            = %e\n', config.MobScrew);
        fprintf(fid, 'MobEdge             = %e\n', config.MobEdge);
        fprintf(fid, 'MobGB               = %f\n\n', config.MobGB);
        fprintf(fid, 'MobClimb            = %e\n', config.MobClimb);
        fprintf(fid, 'MobGlide            = %e\n\n', config.MobGlide);
        fprintf(fid, 'includeInertia      = %d\n', config.includeInertia);
        fprintf(fid, 'massDensity         = %f\n\n', config.massDensity);
    end


    fprintf(fid, '# I/O controls and parameters\n\n');
    fprintf(fid, 'armfile             = %d\n', config.armfile);
    fprintf(fid, 'armfilefreq         = %d\n\n', config.armfilefreq);
    fprintf(fid, 'savecn              = %d\n', config.savecn);
    fprintf(fid, 'savecnfreq          = %d\n', config.savecnfreq);
    fprintf(fid, 'savecncounter       = %d\n\n', config.savecncounter);
    fprintf(fid, 'saveprop            = %d\n', config.saveprop);
    fprintf(fid, 'savepropfreq        = %d\n\n', config.savepropfreq);
    fprintf(fid, 'atomeye             = %d\n', config.atomeye);
    fprintf(fid, 'atomeyefreq         = %d\n', config.atomeyefreq);
    fprintf(fid, 'atomeyesegradius    = %f\n\n', config.atomeyesegradius);
    fprintf(fid, 'writeFlux           = %d\n', config.writeFlux);
    fprintf(fid, 'writeFluxFullDecomp = %d\n', config.writeFluxFullDecomp);
    fprintf(fid, 'writeFluxFullDecompTotals = %d\n', config.writeFluxFullDecompTotals);
    fprintf(fid, 'writeFluxSimpleTotals = %d\n', config.writeFluxSimpleTotals);
    fprintf(fid, 'writeFluxFreq       = %d\n\n', config.writeFluxFreq);


    fprintf(fid, '# Miscellaneous parameters\n\n');
    fprintf(fid, 'FEM_DD_ForceSegTol  = %f\n', config.FEM_DD_ForceSegTol);
    fprintf(fid, 'enforceGlidePlanes  = %d\n', config.enforceGlidePlanes);
    fprintf(fid, 'enableCrossSlip     = %d\n', config.enableCrossSlip);
    fprintf(fid, 'TensionFactor       = %f\n', config.TensionFactor);
    fprintf(fid, 'elasticinteraction  = %d\n', config.elasticinteraction);
    
    fclose(fid);
    fprintf('ParaDiS control file written to: %s\n', fname);
end

function plot_bicrystal(rn, links, box_b, config)
    fprintf('Generating 3D plot of the bicrystal structure...\n');
    figure('Name', 'Generated Bicrystal Structure', 'NumberTitle', 'off', 'Color', 'w');
    ax = axes; hold(ax, 'on');
    gb_x = [-box_b(1)/2, box_b(1)/2, box_b(1)/2, -box_b(1)/2];
    gb_z = [-box_b(3)/2, -box_b(3)/2, box_b(3)/2, box_b(3)/2];
    gb_y = [0, 0, 0, 0];
    patch(ax, gb_x, gb_y, gb_z, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Grain Boundary (y=0)');
    b_vectors = links(:, 3:5); normalized_b = zeros(size(b_vectors));
    for i=1:size(b_vectors,1), vec=b_vectors(i,:); fnz=find(vec,1); if ~isempty(fnz)&&vec(fnz)<0, normalized_b(i,:)=-vec; else, normalized_b(i,:)=vec; end, end
    unique_b_families = unique(round(normalized_b, 8), 'rows');
    num_families = size(unique_b_families, 1); color_palette = lines(num_families);
    legend_handles = gobjects(num_families+1, 1); legend_labels = cell(num_families+1, 1);
    legend_handles(1) = findobj(ax, 'Type', 'Patch'); legend_labels{1} = 'Grain Boundary (y=0)';
    plotted_families = false(num_families, 1);
    for i = 1:size(links, 1)
        p_start=rn(links(i,1),1:3); p_end=rn(links(i,2),1:3);
        family_idx = find(all(abs(unique_b_families - normalized_b(i,:)) < 1e-6, 2));
        if isempty(family_idx), continue; end % safety check
        h = plot3(ax, [p_start(1),p_end(1)], [p_start(2),p_end(2)], [p_start(3),p_end(3)], '-', 'LineWidth', 2.0, 'Color', color_palette(family_idx,:));
        if ~plotted_families(family_idx)
            legend_handles(family_idx+1)=h; legend_labels{family_idx+1}=sprintf('b = %s', mat2str(unique_b_families(family_idx,:), 2));
            plotted_families(family_idx) = true;
        end
    end
    plot3(ax, rn(:,1), rn(:,2), rn(:,3), 'k.', 'MarkerSize', 10);
    axis(ax,'equal'); grid(ax,'on'); box_lims=box_b/2;
    xlim(ax,[-box_lims(1),box_lims(1)]); ylim(ax,[-box_lims(2),box_lims(2)]); zlim(ax,[-box_lims(3),box_lims(3)]);
    xlabel(ax,'X [b]'); ylabel(ax,'Y [b]'); zlabel(ax,'Z [b]');
    title(ax,sprintf('Generated Bicrystal Structure for %s', config.crystalStructure));
    legend_handles = legend_handles(legend_handles~=0);
    legend_labels = legend_labels(~cellfun('isempty',legend_labels));
    legend(ax, legend_handles, legend_labels, 'Location', 'eastoutside', 'FontSize', 9);
    view(ax,3); rotate3d(ax,'on'); hold(ax,'off');
end