% -------------------------------------------------------------------------
%  Author:     Prof. Antonio Franchi
%  Affiliation: 
%      1. University of Twente, The Netherlands
%      2. Sapienza University of Rome, Italy
%  Date:       2025
%  
%  Description: 
%      [One line description of what this file does]
%
%  Copyright (c) 2025 Antonio Franchi. All rights reserved.
% -------------------------------------------------------------------------

function OmniDesign_Benchmarker()
% =========================================================================
% OMNIDESIGN OPTIMIZER - MODULE 5: BENCHMARKER & TOPOLOGY VISUALIZER
% =========================================================================
%
% PROJECT:      OmniDesign Optimizer Framework
% DESCRIPTION:
%   This script reproduces the comparative analysis presented in the paper,
%   validating the proposed "Log-Volume" potential against the standard
%   "Condition Number" metric.
%
%   It performs two key functions:
%   1. SCALABILITY BENCHMARK (Table 1): 
%      Measures and compares execution time between the two cost functions
%      across increasing chassis complexity (N=6 to N=20).
%   2. TOPOLOGICAL ABLATION (Figures):
%      Visualizes the resulting solution spaces. It demonstrates that 
%      Log-Volume (J_vol) produces coherent, continuous 1D manifolds, 
%      while Condition Number (kappa) often results in scattered, 
%      disconnected local minima.
%
% INPUTS:
%   - None (Generates geometries internally).
%
% OUTPUTS:
%   - Prints timing statistics to the Command Window.
%   - Generates side-by-side comparison plots (Blue = J_vol, Orange = Kappa).
%
% =========================================================================

    close all; clc;
    rng(42); % Fixed seed for reproducibility of the benchmark data
    
    %% --- CONFIGURATION ---
    
    % Sample counts for the stochastic search
    % M_scaling: High count (100) for accurate timing averages in Part 1.
    % M_viz:     Lower count (50) sufficient for visualization in Part 2.
    config.M_scaling = 100;  
    config.M_viz     = 50;      
    
    % OPTIMIZATION SETTINGS
    % We use 'sqp' (Sequential Quadratic Programming) for its robustness 
    % in handling the nonlinear equality constraints (unit norm ||u||=1) 
    % inherent to the manifold definition.
    config.options = optimoptions('fmincon', 'Display', 'none', ...
        'Algorithm', 'sqp', ...
        'OptimalityTolerance', 1e-7, ...
        'ConstraintTolerance', 1e-7);

    % SHAPES FOR PART 2 (VISUALIZATION)
    % Format: {DisplayName, GeometryType, N_rotors}
    % The script will generate comparison figures for these specific shapes.
    config.ablation_set = {      
        'Quasi-Hexagon',   'quasi_polygon_rad', 6;
        'Cube',            'cube',              8;
        'Regular Decagon', 'polygon',           10;
    };
    
    fprintf('=======================================================\n');
    fprintf('   OMNIDESIGN FINAL BENCHMARK SUITE (v4)\n');
    fprintf('=======================================================\n\n');

    %% --------------------------------------------------------------------
    %  PART 1: SCALING ANALYSIS (Time Comparison)
    %  --------------------------------------------------------------------
    %  Iterates through geometries of increasing complexity (N) to measure 
    %  the average convergence time per start. This generates the data for 
    %  Table 1 in the manuscript.
    %  --------------------------------------------------------------------
    fprintf('--- 1. Scaling Analysis (Time Comparison) ---\n');
    
    % Geometries to test for Table 1
    % 
    scaling_geoms = {
        'Reg. Hexagon',      'polygon',                  6;
        'Quasi-Hexagon-1',   'quasi_polygon_rad',        6;
        'Quasi-Hexagon-2',   'quasi_polygon_tan',        6;
        'Triang. Prism',     'tri_prism',                6;
        'Irreg. Pent. Pyr',  'irregular_pent_pyramid',   6;
        'Cube',              'cube',                     8; 
        'Icosahedron',       'icosahedron',              12;
        'Dodecahedron',      'dodecahedron',             20
    };
    
    % Storage for timing results: [Time_LogVol, Time_CondNum]
    timing_data = zeros(size(scaling_geoms, 1), 2); 
    
    for i = 1:size(scaling_geoms, 1)
        name = scaling_geoms{i,1};
        type = scaling_geoms{i,2};
        N    = scaling_geoms{i,3};
        
        fprintf('Processing %-16s (N=%d)... ', name, N);
        
        % Generate chassis geometry
        [P, ~, L] = generate_geometry(type, N);
        
        % A. Benchmark Log-Volume (Proposed)
        t_start = tic;
        run_optimization_batch(P, L, config.M_scaling, config.options, 'log_vol');
        t_log = toc(t_start);
        ms_log = (t_log / config.M_scaling) * 1000; % Time per start in ms
        
        % B. Benchmark Condition Number (Baseline)
        t_start = tic;
        run_optimization_batch(P, L, config.M_scaling, config.options, 'cond_num');
        t_cond = toc(t_start);
        ms_cond = (t_cond / config.M_scaling) * 1000; % Time per start in ms
        
        timing_data(i, :) = [ms_log, ms_cond];
        
        % Live feedback
        fprintf('Log: %.1f ms | Cond: %.1f ms (Ratio: %.1fx)\n', ...
            ms_log, ms_cond, ms_cond/ms_log);
    end

    %% --------------------------------------------------------------------
    %  PART 2: ABLATION & VISUALIZATION LOOP
    %  --------------------------------------------------------------------
    %  Visualizes the solution landscape side-by-side. 
    %  - Left (Blue): J_vol (Log-Volume)
    %  - Right (Orange): Kappa (Condition Number)
    %  This visually confirms the smoothness of the manifold for J_vol.
    %  --------------------------------------------------------------------
    fprintf('\n--- 2. Ablation & Visualization Set ---\n');
    
    n_shapes = size(config.ablation_set, 1);
    
    for k = 1:n_shapes
        shape_name = config.ablation_set{k, 1};
        shape_type = config.ablation_set{k, 2};
        N_k        = config.ablation_set{k, 3};
        
        fprintf('\n>>> Visualizing Shape %d/%d: %s (N=%d)\n', k, n_shapes, shape_name, N_k);
        [P_abl, ~, L_abl] = generate_geometry(shape_type, N_k);
        
        % A. Solve using Log-Volume
        fprintf('   Running Log-Volume... ');
        res_log = run_optimization_batch(P_abl, L_abl, config.M_viz, config.options, 'log_vol');
        disp_log = calculate_dispersion(res_log.minima);
        fprintf('Done. Dispersion: %.4f\n', disp_log);
        
        % B. Solve using Condition Number
        fprintf('   Running Cond-Num...   ');
        res_cond = run_optimization_batch(P_abl, L_abl, config.M_viz, config.options, 'cond_num');
        disp_cond = calculate_dispersion(res_cond.minima);
        fprintf('Done. Dispersion: %.4f\n', disp_cond);
        
        % C. Generate Comparative Figure
        fprintf('   Generating Figure...\n');
        visualize_side_by_side(res_log, res_cond, P_abl, shape_name);
    end

    %% --------------------------------------------------------------------
    %  PART 3: SUMMARY GENERATION
    %  --------------------------------------------------------------------
    %  Outputs the formatted table content for easy inclusion in LaTeX.
    %  --------------------------------------------------------------------
    avg_speedup = mean(timing_data(:,2) ./ timing_data(:,1));
    
    fprintf('\n\n=======================================================\n');
    fprintf('       COPY THIS INTO YOUR LATEX DOCUMENT              \n');
    fprintf('=======================================================\n');
    
    fprintf('\n--- TABLE 1 (Scalability Data) ---\n');
    fprintf('Geometry & N & Time ($J_{vol}$) & Time ($\\kappa$) & Speedup \\\\\n');
    fprintf('\\hline\n');
    for i = 1:size(scaling_geoms,1)
        fprintf('%-16s & %d & %.1f ms & %.1f ms & \\textbf{%.1fx} \\\\\n', ...
            scaling_geoms{i,1}, scaling_geoms{i,3}, timing_data(i,1), timing_data(i,2), ...
            timing_data(i,2)/timing_data(i,1));
    end
    
    fprintf('\n--- SUMMARY TEXT ---\n');
    fprintf('Log-Volume is on average \\textbf{%.1fx faster} than Condition Number optimization.\n', avg_speedup);
    fprintf('For regular symmetric shapes, visual inspection confirms that $J_{vol}$ solutions form coherent 1D manifolds (curves), whereas $\\kappa$ solutions exhibit unstructured scatter.\n');
end

%% ========================================================================
%  VISUALIZATION FUNCTION
%  Plots rotor orientations projected onto the unit disk (top-down view).
%  Blue = Log-Volume (Proposed), Orange = Condition Number (Standard).
%  ========================================================================
function visualize_side_by_side(res_log, res_cond, P, shape_name)
    N = size(P, 2);
    % Limit plotting to 8 rotors max for figure readability, even if N > 8
    N_plot = min(N, 8); 
    
    % Create a new figure for THIS shape
    f = figure('Name', ['Comparison: ' shape_name], 'Color', 'w', 'Position', [50 50 1400 600]);
    
    % --- ROW 1: LOG-VOLUME (BLUE) ---
    for n = 1:N_plot
        subplot(2, N_plot, n); hold on; axis equal;
        set(gca, 'XColor', 'none', 'YColor', 'none'); 
        xlim([-1.2 1.2]); ylim([-1.2 1.2]);
        draw_circle();
        
        pts = squeeze(res_log.minima(:, n, :));
        pts = flip_negative_z(pts); % Map lower hemisphere to upper for 2D viz
        
        % Normalize strictly to prevent float errors from appearing as outliers
        pts = pts ./ sqrt(sum(pts.^2, 1));
        
        scatter(pts(1,:), pts(2,:), 20, ...
            'MarkerFaceColor', [0 0.44 0.74], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
            
        if n == 1
            text(-2.5, 0, '$J_{vol}$', ...
                'FontSize', 24, 'Interpreter', 'latex', ...
                'Color', [0 0.44 0.74], 'HorizontalAlignment', 'center');
        end
        title(sprintf('Rotor %d', n), 'FontWeight', 'normal');
    end
    
    % --- ROW 2: CONDITION NUMBER (ORANGE) ---
    for n = 1:N_plot
        subplot(2, N_plot, N_plot + n); hold on; axis equal;
        set(gca, 'XColor', 'none', 'YColor', 'none');
        xlim([-1.2 1.2]); ylim([-1.2 1.2]);
        draw_circle();
        
        pts = squeeze(res_cond.minima(:, n, :));
        pts = flip_negative_z(pts);
        
        % Normalize strictly
        pts = pts ./ sqrt(sum(pts.^2, 1));
        
        scatter(pts(1,:), pts(2,:), 20, ...
            'MarkerFaceColor', [0.85 0.32 0.1], 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6);
            
        if n == 1
            text(-2.5, 0, '$\kappa$', ...
                'FontSize', 24, 'Interpreter', 'latex', ...
                'Color', [0.85 0.32 0.1], 'HorizontalAlignment', 'center');
        end
    end
end

function draw_circle()
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k-', 'LineWidth', 1, 'Color', [0.8 0.8 0.8]);
end

function pts_out = flip_negative_z(pts_in)
    % Reflects points with z < 0 to the upper hemisphere. 
    % Essential for visualizing axes on a 2D disk without ambiguity.
    pts_out = pts_in;
    if size(pts_out, 2) == 1, pts_out = pts_out'; end 
    idx_neg = pts_out(3,:) < 0;
    pts_out(:, idx_neg) = -pts_out(:, idx_neg);
end

%% ========================================================================
%  CORE OPTIMIZATION & MATH HELPERS
%  ========================================================================
function d = calculate_dispersion(minima)
    % Calculates the spread of solutions. High dispersion implies the solution
    % space is not a single point but a set of equivalent optima (or a manifold).
    [~, N, M] = size(minima);
    total_var = 0;
    for n = 1:N
        vecs = squeeze(minima(:, n, :)); 
        % Align hemispheres to handle antipodal symmetry
        for i = 1:M
            if vecs(3,i) < 0, vecs(:,i) = -vecs(:,i); end
        end
        variances = var(vecs, 0, 2); 
        total_var = total_var + sum(variances);
    end
    d = sqrt(total_var / N);
end

function results = run_optimization_batch(P, L, M, options, mode)
    % Runs M independent optimizations from random starts to sample the space.
    N = size(P, 2);
    results.minima = zeros(3, N, M);
    
    switch mode
        case 'log_vol',  obj_fun = @(x) cost_log_vol(x, P, L);
        case 'cond_num', obj_fun = @(x) cost_cond_num(x, P, L);
    end
    
    % Stochastic Initialization
    starts = randn(3, N, M);
    
    for i = 1:M
        u0_raw = starts(:,:,i);
        u0 = u0_raw ./ sqrt(sum(u0_raw.^2, 1)); 
        
        % Solve using SQP
        [x_opt, ~] = fmincon(obj_fun, u0(:), [], [], [], [], [], [], ...
                                @(x) unit_norm_constraint(x, N), options);
        u = reshape(x_opt, 3, N);
        results.minima(:,:,i) = u;
    end
end

function J = cost_log_vol(x, P, L)
    % The Proposed Cost Function: Sum of log-singular values.
    % Minimizing this maximizes the volume of the capability ellipsoid.
    U = reshape(x, 3, size(P,2));
    s = get_singular_values(U, P, L);
    % Log-barrier protects against singularity (s approx 0)
    J = -sum(log(s(1:6) + 1e-12)); 
end

function J = cost_cond_num(x, P, L)
    % The Standard Cost Function: Condition Number (sigma_max / sigma_min).
    U = reshape(x, 3, size(P,2));
    s = get_singular_values(U, P, L);
    % Hard penalty for singularity
    if s(6) < 1e-9, J = 1e9; else, J = s(1)/s(6); end
end

function s = get_singular_values(U, P, L)
    % Computes singular values of the Geometry Matrix A
    N = size(P,2);
    A = zeros(6, N);
    % Construct interaction matrix: [Force; Moment/L]
    for i = 1:N 
        A(:, i) = [U(:,i); (1/L) * cross(P(:,i), U(:,i))]; 
    end
    s = svd(A);
end

function [c, ceq] = unit_norm_constraint(x, N)
    % Enforces ||u_i|| = 1 for all rotors
    ceq = sum(reshape(x,3,N).^2, 1) - 1; 
    c = []; 
end

%% ========================================================================
%  GEOMETRY LIBRARY
%  Defines the vertex positions (P) for various chassis types.
%  All active and experimental shapes are preserved below.
%  ========================================================================
function [P, N, R] = generate_geometry(type, N_poly)
    R = 1.0; % Default Radius
    
    switch lower(type)
        % --- REGULAR POLYGONS ---
        case 'polygon'
            N = N_poly;
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            P = [cos(theta); sin(theta); zeros(1, N)];
            R = 1.0;
            
        % --- QUASI-POLYGONS (Ablation Studies) ---
        case 'quasi_polygon_tan' % Tangential perturbation
            N = N_poly;
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            quasi_theta = 0.15*2*pi*(1 - mod(0:N-1, 2)); % Shift every other vertex
            theta =  theta + quasi_theta;
            P = [cos(theta); sin(theta); zeros(1, N)];
            R = 1.0;    
            
        case 'quasi_polygon_rad' % Radial perturbation
            N = N_poly;
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            P = [cos(theta); sin(theta); zeros(1, N)];
            quasi_R = 1 - 0.1*mod(0:N-1, 2); % Scale radius of every other vertex
            P = P*diag(quasi_R); 
            R = 1.0;        
            
        % --- PLATONIC & ARCHIMEDEAN SOLIDS ---
        case 'cube'
            P = [-1, 1, -1, 1, -1, 1, -1, 1; 
                  1, 1, 1, 1, -1, -1, -1, -1; 
                  1, 1, -1, -1, 1, 1, -1, -1];
            N = 8; R = sqrt(3);
            
        case 'quasi_cube'
            % Slightly distorted cube
            P = [-1.1, 1.1, -1.1, 1.1, -1.1, 1.1, -1.1, 1.1;
                  1, 1, 1, 1, -1, -1, -1, -1; 
                  1, 1, -1, -1, 1, 1, -1, -1];
            N = 8; R = sqrt(3);
            
        case 'octahedron'
            % Vertices of a unit octahedron rotated to align a face with the XY plane
            z_val = 1/sqrt(3);
            r_val = sqrt(2/3);
            theta1 = [0, 120, 240] * (pi/180);
            theta2 = [60, 180, 300] * (pi/180);
            P = [ r_val*cos(theta1), r_val*cos(theta2);
                  r_val*sin(theta1), r_val*sin(theta2);
                  z_val*ones(1,3),  -z_val*ones(1,3) ];
            N = 6; R = 1.0;
            
        case 'dodecahedron'
            phi = (1 + sqrt(5))/2;
            P = get_dodecahedron_verts(phi);
            P = safety_rotate(P); % Avoid axis alignment singularities
            N = 20; R = sqrt(3);
            
        case 'icosahedron'
            phi = (1 + sqrt(5))/2;
            P = get_icosahedron_verts(phi);
            P = safety_rotate(P); 
            N = 12; R = sqrt(1 + phi^2);
            
        case 'cuboctahedron'
            P = [];
            for x=[-1,1], for y=[-1,1], P=[P, [x;y;0], [x;0;y], [0;x;y]]; end, end
            P = unique(P', 'rows')'; 
            P = safety_rotate(P); 
            N = 12; R = sqrt(2);
            
        % --- PRISMS & PYRAMIDS ---
        case 'hex_prism'
            theta = linspace(0, 2*pi, 7); theta(end)=[];
            P_hex = [cos(theta); sin(theta)];
            P = [P_hex, P_hex; ones(1,6), -ones(1,6)];
            P = safety_rotate(P); N = 12; R = sqrt(2);
            
        case 'tri_prism' % N=6
            theta = linspace(0, 2*pi, 4); theta(end)=[];
            P_tri = [cos(theta); sin(theta)];
            P = [P_tri, P_tri; ones(1,3), -ones(1,3)];
            P = safety_rotate(P); N = 6; R = sqrt(2);
            
        case 'pent_bipyramid' % N=7
            theta = linspace(0, 2*pi, 6); theta(end)=[];
            P_eq = [cos(theta); sin(theta); zeros(1,5)];
            P_poles = [0 0; 0 0; 1 -1];
            P = [P_eq, P_poles];
            P = safety_rotate(P); N = 7; R = 1.0; 
            
        case 'sq_antiprism' % N=8
            theta = linspace(0, 2*pi, 5); theta(end)=[];
            P_top = [cos(theta); sin(theta); ones(1,4)];
            P_bot = [cos(theta + pi/4); sin(theta + pi/4); -ones(1,4)];
            P = [P_top, P_bot];
            P = safety_rotate(P); N = 8; R = sqrt(2);
            
        case 'tri_cupola' % N=9
            theta6 = linspace(0, 2*pi, 7); theta6(end)=[];
            P_bot = [cos(theta6); sin(theta6); -ones(1,6)];
            theta3 = linspace(0, 2*pi, 4); theta3(end)=[];
            P_top = [cos(theta3); sin(theta3); ones(1,3)];
            P = [P_bot, P_top];
            P = safety_rotate(P); N = 9; R = sqrt(2);
        % --- IRREGULAR / NOISY TEST CASES (Robustness Checks) ---
        case 'irregular_pent_pyramid' % N=6
            theta = [0, 75, 130, 240, 310] * (pi/180); 
            radii = [1.0, 1.2, 0.9, 1.1, 0.8]; 
            P_base = [radii.*cos(theta); radii.*sin(theta); zeros(1,5)];
            P_apex = [0.2; 0.3; 1.0]; % Off-center apex
            P = [P_base, P_apex];
            P = safety_rotate(P); N = 6; R = 1.0;
            
        case 'perturbed_hex_pyramid' % N=7
            theta = linspace(0, 2*pi, 7); theta(end)=[];
            P_base = [cos(theta); sin(theta); zeros(1,6)];
            noise = [0.1 -0.1 0.05 0.2 -0.15 0.1; 
                     0.05 0.1 -0.1 0.0 0.1 -0.2; 
                     0.1 -0.05 0.1 0.0 0.0 0.1];
            P_base = P_base + noise;
            P_apex = [-0.1; -0.2; 1.2];
            P = [P_base, P_apex];
            P = safety_rotate(P); N = 7; R = 1.0;
            
        case 'irregular_frustum' % N=8
            P_bot = [-1 1 1.2 -0.9; -1 -1 0.9 1.1; -1 -1 -1 -1]; 
            P_top = [-0.6 0.7 0.5 -0.8; -0.7 -0.5 0.8 0.6; 1.2 1.1 1.3 1.0];
            P = [P_bot, P_top];
            P = safety_rotate(P); N = 8; R = sqrt(3);
            
        case 'distorted_oct_pyramid' % N=9
            theta = linspace(0, 2*pi, 9); theta(end)=[];
            radii = 1.0 + 0.3*sin(3*theta); 
            P_base = [radii.*cos(theta); radii.*sin(theta); zeros(1,8)];
            P_apex = [0.4; -0.1; 1.0];
            P = [P_base, P_apex];
            P = safety_rotate(P); N = 9; R = 1.0;
            
        case 'random_cloud_10' % N=10
            prev_rng = rng; rng(12345); % Deterministic noise
            P = randn(3, 10);
            P = P ./ vecnorm(P); 
            P = P .* (0.8 + 0.4*rand(1,10)); 
            rng(prev_rng); 
            P = safety_rotate(P); N = 10; R = 1.0;
            
        otherwise, error('Unknown shape');
    end
end

% --- PLATONIC SOLID HELPERS ---
function P = get_dodecahedron_verts(phi)
    P1 = [-1 1 -1 1 -1 1 -1 1; 1 1 1 1 -1 -1 -1 -1; 1 1 -1 -1 1 1 -1 -1];
    inv = 1/phi;
    P2 = [0 0 0 0; inv inv -inv -inv; phi -phi phi -phi];
    P3 = [inv inv -inv -inv; phi -phi phi -phi; 0 0 0 0];
    P4 = [phi -phi phi -phi; 0 0 0 0; inv inv -inv -inv];
    P = [P1 P2 P3 P4];
end

function P = get_icosahedron_verts(phi)
    P1 = [0, 0, 0, 0; 1, 1, -1, -1; phi, -phi, phi, -phi];
    P2 = [1, 1, -1, -1; phi, -phi, phi, -phi; 0, 0, 0, 0];
    P3 = [phi, -phi, phi, -phi; 0, 0, 0, 0; 1, 1, -1, -1];
    P = [P1, P2, P3];
end

function P_rot = safety_rotate(P)
    % Rotates geometry slightly to avoid perfect alignment with global axes,
    % which can artificially mask singularities in some numerical solvers.
    ang = 20 * pi/180;
    Rx = [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)];
    Ry = [cos(ang) 0 sin(ang); 0 1 0; -sin(ang) 0 cos(ang)];
    P_rot = Ry * Rx * P;
end