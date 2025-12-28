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

function OmniDesign_Optimizer()
% =========================================================================
% OMNIDESIGN OPTIMIZER - MODULE 1: GLOBAL MANIFOLD EXHAUSTION
% =========================================================================
%
% PROJECT:      OmniDesign Optimizer Framework
% DESCRIPTION:  
%   This is the primary entry point for the optimization pipeline. It 
%   performs a "Global Manifold Exhaustion" strategy to find optimal 
%   rotor orientations for a given chassis geometry.
%
%   The algorithm performs the following steps:
%   1. GEOMETRY: Generates the chassis vertex set P (Polygon, Platonic, etc.).
%   2. OPTIMIZATION: Runs a Monte Carlo Multi-Start SQP (Sequential Quadratic 
%      Programming) optimization on the manifold (S^2)^N / Z2.
%   3. PRUNING: Filters results to separate Global Optima (Design Candidates)
%      from Local Minima (Traps).
%   4. SENSITIVITY: (Optional) Sweeps the characteristic length L to 
%      analyze the trade-off between Torque and Force capability.
%   5. EXPORT: Saves the pruned "Global Optima" to a .mat file for 
%      topological analysis in Module 3.
%
% OUTPUT:
%   - .mat file containing the 'PointsTable' (optimal configurations).
%   - Visualization of the solution landscape (RP2 discs).
%
% =========================================================================

    close all; 
    clc;
    
    %% ====================================================================
    %  1. USER CONFIGURATION
    %  ====================================================================
    
    % --- Chassis Selection ---
    % Defines the spatial arrangement of the rotor centers (chassis geometry).
    % Options:
    %   Planar:   'polygon', 'quasi_polygon_rad', 'quasi_polygon_tan'
    %   Platonic: 'octahedron', 'cube', 'dodecahedron', 'icosahedron'
    %   Prisms:   'tri_prism', 'hex_prism', 'sq_antiprism'
    %   Irregular: 'pent_bipyramid', 'tri_cupola', 'quasi_cube', etc.
    config.shape = 'polygon'; 
    
    % --- Optimization Parameters ---
    config.polygon_N = 9;       % Number of rotors (Only used if shape is 'polygon')
    config.M         = 1000;    % Monte Carlo sample size (Higher M = denser landscape map)
    config.L_nominal = 'auto';  % Characteristic Length scale ('auto' = bounding radius)
    
    % --- Pruning Parameters ---
    % The tolerance within which a solution is considered a "Global Minimum"
    % relative to the best found cost.
    config.pruning_tol = 1e-4;  
    
    % --- Sensitivity Analysis ---
    % If true, performs a sweep of L to check conditioning stability.
    config.run_sensitivity = true;
    config.L_range = logspace(-1, 1, 20); 
    
    % --- Visualization Settings ---
    cfg.geom_scale = 1.5; 
    cfg.font_scale = 1.4;
    
    % --- Solver Options (SQP) ---
    % Using 'fmincon' with SQP is efficient for the smooth constraints of S^2.
    config.options = optimoptions('fmincon', ...
        'Display', 'none', ...
        'Algorithm', 'sqp', ...
        'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-7);

    %% ====================================================================
    %  2. GEOMETRY GENERATION
    %  ====================================================================
    fprintf('=== STEP 1: Geometry Generation ===\n');
    
    % Generate the chassis vertex set P (3xN matrix)
    [P, N_propellers, shape_radius] = generate_geometry(config);
    config.N = N_propellers;
    
    % Set the characteristic length if 'auto'
    if strcmp(config.L_nominal, 'auto')
        config.L_nominal = shape_radius;
    end
    fprintf('Shape: %s (N=%d) | L = %.4f\n', config.shape, config.N, config.L_nominal);
    
    % Debug Plot: Visual check of the chassis geometry
    plotPoly3(P);

    %% ====================================================================
    %  3. GLOBAL SEARCH & TOPOLOGICAL PRUNING
    %  ====================================================================
    fprintf('\n=== STEP 2: Main Optimization (M=%d samples) ===\n', config.M);
    
    % 1. Run Global Search (Monte Carlo)
    %    Explores the landscape from M random starting points.
    results_main = run_optimization_batch(P, config.L_nominal, config.M, config.options);
    
    % 2. PRUNING LOGIC
    %    We separate the "True Solutions" (Global Minima) from "Traps" (Local Minima).
    min_cost = min(results_main.costs);
    
    % Indices of Global Optima (within numerical tolerance)
    idx_global = find(results_main.costs <= (min_cost + config.pruning_tol));
    
    % Indices of Local Minima (suboptimal traps)
    idx_local  = setdiff(1:config.M, idx_global);
    
    % --- DETAILED PRUNING REPORT ---
    fprintf('\n------------------------------------------------------------\n');
    fprintf('                 GLOBAL SEARCH PRUNING REPORT               \n');
    fprintf('------------------------------------------------------------\n');
    fprintf('Total Monte Carlo Starts  : %d\n', config.M);
    fprintf('Global Optima Retained    : %d (%.1f%%)\n', length(idx_global), 100*length(idx_global)/config.M);
    fprintf('Local Minima Discarded    : %d (%.1f%%)\n', length(idx_local), 100*length(idx_local)/config.M);
    fprintf('------------------------------------------------------------\n');
    fprintf('Global Minimum Cost       : %.8f\n', min_cost);
    
    if ~isempty(idx_local)
        costs_local = results_main.costs(idx_local);
        best_local = min(costs_local);
        worst_local = max(costs_local);
        avg_local = mean(costs_local);
        
        % The "Gap" quantifies the energy barrier between the optimal manifold
        % and the nearest local trap.
        fprintf('Best Local Minimum Cost   : %.8f (Gap: %.2e)\n', best_local, best_local - min_cost);
        fprintf('Worst Local Minimum Cost  : %.8f\n', worst_local);
        fprintf('Avg Local Cost            : %.8f\n', avg_local);
    else
        fprintf('Convergence               : PERFECT (No local minima found)\n');
    end
    fprintf('------------------------------------------------------------\n');

    %% ====================================================================
    %  4. SENSITIVITY ANALYSIS (Optional)
    %  ====================================================================
    results_sens = [];
    if config.run_sensitivity
        fprintf('\n=== STEP 3: Sensitivity Analysis ===\n');
        % Sweeps L to observe the Condition Number (Isotropy) vs. Singular Value (Strength)
        results_sens = run_sensitivity_sweep(P, config.L_range, config.options);
    end

    %% ====================================================================
    %  5. DATA EXPORT
    %  ====================================================================
    % Save the results. 
    % NOTE: 'PointsTable' will contain ONLY the pruned Global Optima.
    % This ensures the subsequent Topological Analysis module receives clean data.
    save_data_2D_pruned(results_main, idx_global, results_sens, P, config);

    %% ====================================================================
    %  6. VISUALIZATION (DUAL LAYER)
    %  ====================================================================
    fprintf('\n=== STEP 4: Generating Figures ===\n');
    
    % Viz 1: The Initial Random Distribution (Uniformity Check)
    visualize_starts(results_main, P, 'Initial Starting Points', ...
                     cfg.geom_scale, cfg.font_scale);
                     
    % Viz 2: The Solution Landscape (RP2)
    %        Displays Global Optima (Solid) overlaid on Local Minima (Faint/Ghost).
    visualize_solutions_layered(results_main, idx_global, idx_local, P, ...
                        'Solution Landscape (RP2)', cfg.geom_scale, cfg.font_scale);
    
    % Viz 3: Sensitivity Plots (if enabled)
    if config.run_sensitivity
        visualize_sensitivity(results_sens, config.L_nominal, ...
                              'Sensitivity Analysis', cfg.geom_scale, cfg.font_scale);
    end
    
    fprintf('\nDone.\n');
end

%% ========================================================================
%  CORE ALGORITHMS & SOLVERS
%  ========================================================================

function results = run_optimization_batch(P, L, M, options)
% Runs the SQP optimizer M times from random starting points.
    N = size(P, 2);
    results.starts = zeros(3, N, M);
    results.minima = zeros(3, N, M);
    results.costs  = zeros(M, 1);
    results.cond   = zeros(M, 1);
    
    fprintf('Progress: ');
    for i = 1:M
        if mod(i, 50) == 0, fprintf('.'); end
        if mod(i, 500) == 0, fprintf('%d', i); end
        
        % 1. Stochastic Start (Random point on S^2 sphere)
        u0_raw = randn(3, N);
        u0 = u0_raw ./ sqrt(sum(u0_raw.^2, 1)); 
        results.starts(:,:,i) = u0;
        
        % 2. SQP Optimization
        objective = @(x) cost_function_log_vol(x, P, L);
        
        % fmincon with non-linear equality constraint (unit norm)
        [x_opt, fval] = fmincon(objective, u0(:), [], [], [], [], [], [], ...
                                @(x) unit_norm_constraint(x, N), options);
        
        u_opt = reshape(x_opt, 3, N);
        results.minima(:,:,i) = u_opt;
        results.costs(i) = fval;
        
        % 3. Calculate Condition Number of the result
        s_final = get_singular_values(u_opt, P, L);
        results.cond(i) = s_final(1) / s_final(6);
    end
    fprintf(' Done.\n');
end

function sens_data = run_sensitivity_sweep(P, L_range, options)
% Sweeps the characteristic length L to analyze matrix properties.
    N = size(P, 2);
    n_steps = length(L_range);
    sens_data.L = L_range;
    sens_data.cond = zeros(n_steps, 1);
    sens_data.sigma_min = zeros(n_steps, 1);
    
    fprintf('Sweeping L: ');
    for k = 1:n_steps
        L_curr = L_range(k);
        if mod(k, 5) == 0, fprintf('.'); end
        
        best_fval = inf; 
        best_x = [];
        % Perform a mini multi-start (3 attempts) to ensure we find the 
        % global optimum for this specific L value.
        for attempt = 1:3 
            u0 = randn(3, N); u0 = u0./vecnorm(u0);
            [x_opt, fval] = fmincon(@(x) cost_function_log_vol(x, P, L_curr), ...
                u0(:), [], [], [], [], [], [], ...
                @(x) unit_norm_constraint(x, N), options);
            if fval < best_fval, best_fval = fval; best_x = x_opt; end
        end
        U = reshape(best_x, 3, N);
        s = get_singular_values(U, P, L_curr);
        
        sens_data.cond(k) = s(1)/s(6);      % Condition Number (Isotropy)
        sens_data.sigma_min(k) = s(6);      % Minimum Singular Value (Worst-case Authority)
    end
    fprintf(' Done.\n');
end

%% ========================================================================
%  COST FUNCTION & CONSTRAINTS
%  =======================================================================

function J = cost_function_log_vol(x, P, L)
% Cost: Negative Log Volume of the wrench polytope.
% J = -sum( log( singular_values ) )
% This is equivalent to maximizing the volume of the ellipsoid.
    U = reshape(x, 3, size(P,2));
    s = get_singular_values(U, P, L);
    % Add epsilon to prevent log(0)
    J = -sum(log(s(1:6) + 1e-12)); 
end

function s = get_singular_values(U, P, L)
% Computes singular values of the Geometry Matrix A.
    N = size(P,2);
    A = zeros(6, N);
    for i = 1:N
        % Column i: [Direction; Moment/L]
        A(:, i) = [U(:,i); (1/L) * cross(P(:,i), U(:,i))]; 
    end
    s = svd(A);
end

function [c, ceq] = unit_norm_constraint(x, N)
% Non-linear equality constraint: ||u_i|| = 1 for all i
    ceq = sum(reshape(x,3,N).^2, 1) - 1; 
    c = []; 
end

%% ========================================================================
%  VISUALIZATION FUNCTIONS
%  ========================================================================

function visualize_solutions_layered(results, idx_global, idx_local, P, title_str, geom_scale, font_scale)
% Plots solutions on the unit disc (RP2 representation).
% LAYER 1 (Background): Local Minima (Faint, "Ghosts")
% LAYER 2 (Foreground): Global Optima (Solid, "Real")
    
    N = size(P, 2);
    fs = 14 * font_scale;
    
    % --- Visualization Style ---
    col_local  = [0.6 0.6 0.7]; % Faint Blue-Grey
    alpha_local = 0.15;         % Very transparent
    
    col_global = [0.85 0.32 0.1]; % Highlight Color (Burnt Orange/Red)
    alpha_global = 0.8;           % Solid
    
    f = figure('Name', title_str, 'Color', 'w', 'Position', [100 100 1200 700]);
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'tight');
    title(t, sprintf('\\textbf{%s}', title_str), 'Interpreter', 'latex', 'FontSize', fs + 4);
    
    for n = 1:N
        nexttile; hold on; axis equal; 
        xlim([-1.1, 1.1]); ylim([-1.1, 1.1]);
        box off; axis off;
        
        % Draw Disc Boundary (Unit Circle)
        theta = linspace(0, 2*pi, 200);
        plot(cos(theta), sin(theta), 'k-', 'LineWidth', 0.8 * geom_scale);
        % Draw Crosshairs
        plot([-1 1], [0 0], ':', 'Color', [0.9 0.9 0.9]);
        plot([0 0], [-1 1], ':', 'Color', [0.9 0.9 0.9]);
        
        % --- PROCESS POINTS (Canonical Mapping z > 0) ---
        pts_all = results.minima(:,:,:);
        
        % Extract mapped points (x,y)
        [x_loc, y_loc] = extract_and_map(pts_all, n, idx_local);
        [x_glob, y_glob] = extract_and_map(pts_all, n, idx_global);
        
        % --- LAYER 1: LOCAL MINIMA ---
        if ~isempty(x_loc)
            scatter(x_loc, y_loc, 15 * geom_scale, ...
                'MarkerFaceColor', col_local, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', alpha_local);
        end
        
        % --- LAYER 2: GLOBAL OPTIMA ---
        if ~isempty(x_glob)
            scatter(x_glob, y_glob, 20 * geom_scale, ...
                'MarkerFaceColor', col_global, ...
                'MarkerEdgeColor', 'k', ... 
                'LineWidth', 0.5, ...
                'MarkerFaceAlpha', alpha_global);
        end
        
        % Label
        text(0.01, 0.999, sprintf('Rotor %d', n), 'Units', 'normalized', ...
            'FontSize', fs, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end
    
    % Legend
    hL = scatter(NaN,NaN, 'MarkerFaceColor', col_local, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
    hG = scatter(NaN,NaN, 'MarkerFaceColor', col_global, 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', 'k');
    legend([hL, hG], {'Local Minima', 'Global Optima'}, 'Location', 'bestoutside', 'Orientation','horizontal');
end

function [x_out, y_out] = extract_and_map(pts_all, rotor_idx, sample_indices)
% Helper: Extracts vectors and performs RP2 mapping (flips if z < 0).
    if isempty(sample_indices)
        x_out = []; y_out = []; return;
    end
    
    vecs = squeeze(pts_all(:, rotor_idx, sample_indices)); % 3 x M_subset
    
    % Handle dimensions for single-point cases
    if size(vecs, 2) == 1 && size(vecs, 1) == 3
        % Correct orientation
    elseif size(vecs, 1) ~= 3
        vecs = vecs'; 
    end
    
    % Canonical Mapping: If z < 0, flip vector to the upper hemisphere
    neg_z_mask = vecs(3,:) < 0;
    vecs(:, neg_z_mask) = -vecs(:, neg_z_mask);
    
    x_out = vecs(1, :);
    y_out = vecs(2, :);
end

function visualize_starts(results, P, title_str, geom_scale, font_scale)
% Plots the initial random starting points to verify uniform coverage.
    N = size(P, 2);
    fs = 14 * font_scale;
    f = figure('Name', title_str, 'Color', 'w', 'Position', [100 100 1000 600]);
    t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'tight');
    title(t, sprintf('\\textbf{%s}', title_str), 'Interpreter', 'latex', 'FontSize', fs+2);
    
    for n = 1:N
        nexttile; hold on; axis equal; axis off;
        theta = linspace(0, 2*pi, 200);
        plot(cos(theta), sin(theta), 'k-', 'LineWidth', 0.5);
        
        pts = squeeze(results.starts(:, n, :));
        scatter(pts(1,:), pts(2,:), 10 * geom_scale, ...
                'MarkerFaceColor', [0 0.4 0.6], 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 0.2);
        title(sprintf('Rotor %d', n), 'FontSize', fs);
    end
end

function visualize_sensitivity(sens, L_nom, title_str, geom_scale, font_scale)
% Plots Condition Number and Min Singular Value vs Characteristic Length L.
    fs = 16 * font_scale;
    figure('Name', title_str, 'Color', 'w', 'Position', [100 100 900 400]);
    
    subplot(1, 2, 1); hold on; grid on;
    loglog(sens.L, sens.cond, 'o-', 'LineWidth', 1.5, 'Color', [0 0.44 0.74]);
    xline(L_nom, '--r');
    xlabel('L'); ylabel('\kappa'); title('Isotropy'); set(gca, 'FontSize', fs);
    
    subplot(1, 2, 2); hold on; grid on;
    semilogx(sens.L, sens.sigma_min, 'o-', 'LineWidth', 1.5, 'Color', [0.85 0.32 0.1]);
    xline(L_nom, '--r');
    xlabel('L'); ylabel('\sigma_{min}'); title('Strength'); set(gca, 'FontSize', fs);
end

%% ========================================================================
%  DATA EXPORT (PRUNED ONLY)
%  ========================================================================

function save_data_2D_pruned(res_main, idx_global, res_sens, P, cfg)
% Exports results to a .mat file.
% Creates a table 'PointsTable' containing ONLY Global Optima.
    N = size(P, 2);
    num_global = length(idx_global);
    
    % Prepare table for Global Optima only
    raw_results = zeros(num_global, 2*N + 2);
    
    for k = 1:num_global
        orig_idx = idx_global(k);
        u_mat = res_main.minima(:,:,orig_idx);
        
        % Canonical Mapping ensures all vectors are in Upper Hemisphere
        for r=1:N
            if u_mat(3,r) < 0, u_mat(:,r) = -u_mat(:,r); end
        end
        
        xy_pairs = u_mat(1:2, :);
        raw_results(k, 1:2*N) = xy_pairs(:)';
        
        % We save the negative condition number and the raw cost
        raw_results(k, end-1) = -res_main.cond(orig_idx); 
        raw_results(k, end)   = res_main.costs(orig_idx);
    end
    
    varNames = {};
    for k = 1:N, varNames{end+1} = sprintf('Prop%d_x', k); varNames{end+1} = sprintf('Prop%d_y', k); end
    varNames{end+1} = 'MinusCond'; varNames{end+1} = 'Cost';
    
    PointsTable = array2table(raw_results, 'VariableNames', varNames);
    filename = sprintf('Results_%s_N%d_GlobalPruned.mat', cfg.shape, N);
    
    % Save:
    % 1. PointsTable: Clean data for analysis
    % 2. res_main: Full raw data (including local minima) for audit/debug
    % 3. P, cfg: Metadata
    save(filename, 'PointsTable', 'res_main', 'res_sens', 'P', 'cfg');
    fprintf('Results saved to %s (Table contains only Global Optima).\n', filename);
end

%% ========================================================================
%  GEOMETRY DEFINITIONS
%  ========================================================================

function [P, N, R] = generate_geometry(cfg)
% GENERATE_GEOMETRY
% Constructs the vertex positions P for the requested chassis shape.
% Returns:
%   P: 3xN matrix of rotor positions
%   N: Number of rotors
%   R: Bounding radius (used as default L)

    switch lower(cfg.shape)
        
        % --- 1. REGULAR POLYGONS (PLANAR) ---
        case 'polygon'
            N = cfg.polygon_N;
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            P = [cos(theta); sin(theta); zeros(1, N)];
            R = 1.0;
            
        case 'quasi_polygon_tan'
            % Polygon with Tangential Perturbation (breaks symmetry)
            N = cfg.polygon_N;
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            quasi_theta = 0.15*2*pi*(1 - mod(0:N-1, 2)); % Zig-zag shift
            theta =  theta + quasi_theta;
            P = [cos(theta); sin(theta); zeros(1, N)];
            R = 1.0;    
            
        case 'quasi_polygon_rad'
            % Polygon with Radial Perturbation (breaks symmetry)
            N = cfg.polygon_N;
            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            P = [cos(theta); sin(theta); zeros(1, N)];
            quasi_R = 1 - 0.1*mod(0:N-1, 2); % In-out radius
            P = P*diag(quasi_R); 
            R = 1.0;        

        % --- 2. PLATONIC SOLIDS ---
        case 'cube'
            % Standard Cube
            P = [-1, 1, -1, 1, -1, 1, -1, 1; 
                  1, 1, 1, 1, -1, -1, -1, -1; 
                  1, 1, -1, -1, 1, 1, -1, -1];
            N = 8; R = sqrt(3);

        case 'quasi_cube'
            % Cube with distortion to break perfect symmetry
            P = [-1.1, 1.1, -1.1, 1.1, -1.1, 1.1, -1.1, 1.1;
                  1, 1, 1, 1, -1, -1, -1, -1; 
                  1, 1, -1, -1, 1, 1, -1, -1];
            N = 8; R = sqrt(3);

        case 'octahedron'
            % Unit octahedron
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
            P = safety_rotate(P); 
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

        % --- 3. PRISMS AND PYRAMIDS ---
        case 'hex_prism'
            theta = linspace(0, 2*pi, 7); theta(end)=[];
            P_hex = [cos(theta); sin(theta)];
            P = [P_hex, P_hex; ones(1,6), -ones(1,6)];
            P = safety_rotate(P); 
            N = 12; R = sqrt(2);
            
        case 'tri_prism' 
            theta = linspace(0, 2*pi, 4); theta(end)=[];
            P_tri = [cos(theta); sin(theta)];
            P = [P_tri, P_tri; ones(1,3), -ones(1,3)];
            P = safety_rotate(P); 
            N = 6; R = sqrt(2);
            
        case 'pent_bipyramid' 
            theta = linspace(0, 2*pi, 6); theta(end)=[];
            P_eq = [cos(theta); sin(theta); zeros(1,5)];
            P_poles = [0 0; 0 0; 1 -1];
            P = [P_eq, P_poles];
            P = safety_rotate(P); 
            N = 7; R = 1.0; 
            
        case 'sq_antiprism' 
            theta = linspace(0, 2*pi, 5); theta(end)=[];
            P_top = [cos(theta); sin(theta); ones(1,4)];
            P_bot = [cos(theta + pi/4); sin(theta + pi/4); -ones(1,4)];
            P = [P_top, P_bot];
            P = safety_rotate(P); 
            N = 8; R = sqrt(2);
            
        case 'tri_cupola' 
            theta6 = linspace(0, 2*pi, 7); theta6(end)=[];
            P_bot = [cos(theta6); sin(theta6); -ones(1,6)];
            theta3 = linspace(0, 2*pi, 4); theta3(end)=[];
            P_top = [cos(theta3); sin(theta3); ones(1,3)];
            P = [P_bot, P_top];
            P = safety_rotate(P); 
            N = 9; R = sqrt(2);
            
        % --- 4. IRREGULAR TEST CASES (Noise/Distortion) ---
        case 'irregular_pent_pyramid' 
            theta = [0, 75, 130, 240, 310] * (pi/180); 
            radii = [1.0, 1.2, 0.9, 1.1, 0.8]; 
            P_base = [radii.*cos(theta); radii.*sin(theta); zeros(1,5)];
            P_apex = [0.2; 0.3; 1.0]; 
            P = [P_base, P_apex];
            P = safety_rotate(P); N = 6; R = 1.0;
            
        case 'perturbed_hex_pyramid' 
            theta = linspace(0, 2*pi, 7); theta(end)=[];
            P_base = [cos(theta); sin(theta); zeros(1,6)];
            noise = [0.1 -0.1 0.05 0.2 -0.15 0.1; 
                     0.05 0.1 -0.1 0.0 0.1 -0.2; 
                     0.1 -0.05 0.1 0.0 0.0 0.1];
            P_base = P_base + noise;
            P_apex = [-0.1; -0.2; 1.2];
            P = [P_base, P_apex];
            P = safety_rotate(P); N = 7; R = 1.0;
            
        case 'irregular_frustum' 
            P_bot = [-1 1 1.2 -0.9; -1 -1 0.9 1.1; -1 -1 -1 -1]; 
            P_top = [-0.6 0.7 0.5 -0.8; -0.7 -0.5 0.8 0.6; 1.2 1.1 1.3 1.0];
            P = [P_bot, P_top];
            P = safety_rotate(P); N = 8; R = sqrt(3);
            
        case 'distorted_oct_pyramid' 
            theta = linspace(0, 2*pi, 9); theta(end)=[];
            radii = 1.0 + 0.3*sin(3*theta); 
            P_base = [radii.*cos(theta); radii.*sin(theta); zeros(1,8)];
            P_apex = [0.4; -0.1; 1.0];
            P = [P_base, P_apex];
            P = safety_rotate(P); N = 9; R = 1.0;
            
        case 'random_cloud_10' 
            prev_rng = rng; rng(12345); % Deterministic Seed
            P = randn(3, 10);
            P = P ./ vecnorm(P); 
            P = P .* (0.8 + 0.4*rand(1,10)); 
            rng(prev_rng); 
            P = safety_rotate(P); N = 10; R = 1.0;
            
        otherwise, error('Unknown shape');
    end
end

function P = get_dodecahedron_verts(phi)
% Helper to generate standard dodecahedron vertices
    P1 = [-1, 1, -1, 1, -1, 1, -1, 1; 1, 1, 1, 1, -1, -1, -1, -1; 1, 1, -1, -1, 1, 1, -1, -1];
    inv_phi = 1/phi;
    P2 = [0, 0, 0, 0; inv_phi, inv_phi, -inv_phi, -inv_phi; phi, -phi, phi, -phi];
    P3 = [inv_phi, inv_phi, -inv_phi, -inv_phi; phi, -phi, phi, -phi; 0, 0, 0, 0];
    P4 = [phi, -phi, phi, -phi; 0, 0, 0, 0; inv_phi, inv_phi, -inv_phi, -inv_phi];
    P = [P1, P2, P3, P4];
end

function P_rot = safety_rotate(P)
% Applies a small rotation to chassis coordinates to avoid 
% numerical alignment with principal axes during optimization.
    ang = 20 * pi/180;
    Rx = [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)];
    Ry = [cos(ang) 0 sin(ang); 0 1 0; -sin(ang) 0 cos(ang)];
    P_rot = Ry * Rx * P;
end

function plotPoly3(P)
% Simple 3D visualization of the chassis geometry
    figure; hold on; axis equal; grid on; view(3);
    
    % Check if the shape is effectively 2D (coplanar on Z)
    if range(P(3,:)) < 1e-6
        % --- 2D CASE (Planar Polygon) ---
        if size(P,2) >= 3
            try
                % Calculate 2D hull using only X and Y
                k = convhull(P(1,:), P(2,:)); 
                % Plot as a flat patch at the Z-level
                patch(P(1,k), P(2,k), P(3,k), 'c', 'FaceAlpha', 0.3);
            catch
                % Fallback for collinear points
            end
        end
    else
        % --- 3D CASE (Polyhedron) ---
        if size(P,2) >= 4
            try
                K = convhull(P(1,:), P(2,:), P(3,:)); 
                trisurf(K, P(1,:), P(2,:), P(3,:), ...
                    'FaceColor','c', 'FaceAlpha',0.3, 'EdgeColor','none'); 
            catch
                warning('Points are coplanar or collinear in 3D; skipping hull.');
            end
        end
    end
    
    % Draw the vertices and edges
    plot3(P(1,:), P(2,:), P(3,:), 'ko', 'MarkerFaceColor','k');
    
    % For polygons, connect the loop
    if range(P(3,:)) < 1e-6
        plot3([P(1,:) P(1,1)], [P(2,:) P(2,1)], [P(3,:) P(3,1)], 'k-');
    end
end