function OmniDesign_Analysis_N5Law_Nfrac_approx()
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

% =========================================================================
% OMNIDESIGN OPTIMIZER - MODULE 3: TOPOLOGICAL ANALYSIS (THE N-5 LAW)
% =========================================================================
%
% PROJECT:      OmniDesign Optimizer Framework
% DESCRIPTION:
%   This is the core analytical engine of the framework. It processes the 
%   intrinsic coordinates extracted in Module 2 to uncover the fundamental 
%   topological structure of the solution space (the "Manifold").
%
%   Key Theoretical Concepts:
%   1. CHIRALITY (Slope): Determines if a rotor rotates CW or CCW relative 
%      to a master reference in the solution space.
%   2. RELATIVE PHASE SPACE: By subtracting the linear motion (Slope * Theta_ref),
%      we collapse the 1-dimensional solution curves into 0-dimensional 
%      points (Clusters). Each cluster represents a distinct physical 
%      Isomer (Topological Branch).
%   3. THE N-5 LAW: Validates the prediction that for N rotors, there are 
%      typically N-5 stable isomers.
%   4. FRACTIONAL LOCKING: Checks if the phase shifts lock to rational 
%      fractions of Pi (specifically multiples of pi/N).
%
% METHODOLOGY:
%   1. Loads 'TorusMapping_Data.mat' (Intrinsic Thetas).
%   2. Performs PCA to determine local tangent slopes (Chirality).
%   3. Transforms data to Relative Phase Space.
%   4. Performs Unsupervised Clustering (K-Means + Silhouette) to identify Isomers.
%   5. Extracts precise Affine Parameters (Phase Shift, Tightness/Spread).
%   6. Fits phases to theoretical N-fractions (k/N * pi).
%   7. Visualizes the Manifold on 2D Flat Torus (Square) or 3D Torus.
%
% INPUT:
%   - 'TorusMapping_Data.mat' (from Module 2).
%
% OUTPUT:
%   - Appends 'IsomerParameters' and 'T_Isomers' table to the .mat file.
%   - Generates visualization figures.
% =========================================================================

    %% 1. Load Pre-Processed Data
    %  Loads the intrinsic angular coordinates computed in Module 2.
    if exist('TorusMapping_Data.mat', 'file')
        dataFile = 'TorusMapping_Data.mat';
    else
        [file, path] = uigetfile('*.mat', 'Select TorusMapping_Data.mat');
        if isequal(file,0), return; end
        dataFile = fullfile(path, file);
    end
    load(dataFile);
    fprintf('Loaded Torus Data for N=%d Disks (%d solutions).\n', N, n_kept);

    %% 2. Configuration UI
    %  Allows the user to select which rotors to visualize as "Reference Axes"
    %  and choose between 2D (Flat Torus) or 3D (Donut) visualization modes.
    M = ceil(sqrt(N)); 
    fig_height = 200 + (M * 30);
    d = dialog('Position',[300 300 400 fig_height], 'Name', 'Step 2: N-5 Analysis Config');
    
    uicontrol('Parent',d,'Style','text','Position',[20 fig_height-30 350 20],...
        'String','Select Master Disks for Correlation Plots:','HorizontalAlignment','left');
        
    chk_disks = gobjects(1,N);
    grid_start_y = fig_height - 60;
    
    % Checkbox Grid Generation
    for i = 1:N
        row = ceil(i/M); col = mod(i-1, M) + 1;
        x_pos = 20 + (col-1)*60;
        y_pos = grid_start_y - (row-1)*30;
        chk_disks(i) = uicontrol('Parent',d,'Style','checkbox',...
            'Position',[x_pos y_pos 50 20],...
            'String',sprintf('D%d',i), 'Value', (i==1));
    end
    
    uicontrol('Parent',d,'Style','text','Position',[20 80 350 20],...
        'String','Visualization Mode:','HorizontalAlignment','left');
    bg = uibuttongroup('Parent',d,'Position',[0.05 0.15 0.9 0.12],'BorderType','none');
    rb_2d = uicontrol('Parent',bg,'Style','radiobutton','String','2D Linear (Square)',...
        'Position',[10 5 150 20]);
    rb_3d = uicontrol('Parent',bg,'Style','radiobutton','String','3D Donut (Torus)',...
        'Position',[170 5 150 20]);
        
    uicontrol('Parent',d,'Position',[150 10 100 30],'String','Analyze',...
        'Callback','uiresume(gcbf)');
    uiwait(d);
    
    if ~ishandle(d), return; end
    
    % Retrieve Selection
    selectedMasters = [];
    for k = 1:N
        if get(chk_disks(k),'Value'), selectedMasters(end+1) = k; end %#ok<AGROW>
    end
    use3D = get(rb_3d, 'Value');
    delete(d);

    %% 3. Global ND Structural Fitting (Local Tangent Method)
    %  This section reduces the manifold topology.
    fprintf('\n--- STARTING TOPOLOGICAL ANALYSIS ---\n');
    
    % A. Determine Slopes (Chirality)
    %    Using PCA on local neighborhoods to find the tangent vector direction.
    masterRef = 1; % Usually Rotor 1 is the reference
    Slopes = determineSlopes(AngularCoords, masterRef, N);
    
    % B. Transform to Relative Phase Space
    %    Subtracting the linear correlation (Slope * MasterPhase) leaves only
    %    the constant phase offset. This collapses the lines into points (Clusters).
    RelativePhases = zeros(n_kept, N);
    for k = 1:N
        diffs = AngularCoords(:, k) - Slopes(k) * AngularCoords(:, masterRef);
        RelativePhases(:, k) = mod(diffs, pi);
    end
    
    % C. Clustering (Detecting Isomers)
    %    We use unsupervised K-means with Silhouette analysis to automatically
    %    determine how many distinct branches (Isomers) exist.
    [bestK, ClusterIdx, ClusterCentroids] = performClustering(RelativePhases, masterRef, n_kept, N);
    
    fprintf('Topological Result: Identified %d Disconnected Branches (Isomers).\n', bestK);
    expectedK = max(1, N-5);
    if N >= 6 && N <= 10
        fprintf('Theoretical Prediction (N-5 Law): K = %d. Match? %d\n', expectedK, (bestK == expectedK));
    end

    %% 4. PARAMETER EXTRACTION & REPORTING (TABLE GENERATION)
    %  Calculates detailed statistics for each identified Isomer, including
    %  physical tightness (Spread) and Fractional Phase Locking.
    fprintf('\n======================================================================================\n');
    fprintf('   EXTRACTED AFFINE PARAMETERS & FRACTIONAL FIT CONFIDENCE  \n');
    fprintf('======================================================================================\n');
    
    % Initialize storage structure
    IsomerParameters = struct('BranchID', {}, 'RotorID', {}, 'Chirality', {}, ...
                              'PhaseShift_Rad', {}, 'PhaseShift_Deg', {}, ...
                              'ClusterRMS_Deg', {}, 'N_Fraction_Str', {}, 'Fit_Error_Deg', {});
    
    % Prepare arrays for a clean Table object
    Tab_Branch = []; Tab_Rotor = []; Tab_Slope = []; Tab_Deg = []; 
    Tab_Spread = []; Tab_NFrac = {}; Tab_FitErr = [];
    
    for b = 1:bestK
        % --- Calculate Global Isomer Confidence (Spread) ---
        % Measures how "tight" the physical solution is around the mathematical center.
        branch_mask = (ClusterIdx == b);
        cluster_points = RelativePhases(branch_mask, :); 
        centroid = ClusterCentroids(b, :); 
        
        % N-dimensional Toroidal Distance RMS (handling the cyclic boundary pi)
        diffs = abs(cluster_points - centroid);
        diffs = min(diffs, pi - diffs); 
        sq_dists = sum(diffs.^2, 2); 
        iso_spread_deg = rad2deg(sqrt(mean(sq_dists)));
        
        fprintf('\n--- Branch (Isomer) #%d [Tightness: +/- %.2f deg] ---\n', b, iso_spread_deg);
        % Table Header
        fprintf('| Rotor | Slope | Phase (Deg) | IsoSpread | Approx (Frac) | Fit Err (deg) |\n');
        fprintf('|-------|-------|-------------|-----------|---------------|---------------|\n');
        
        for k = 1:N
            raw_phase = ClusterCentroids(b, k);
            deg_val = rad2deg(raw_phase);
            slope_val = Slopes(k);
            
            % --- FRACTION APPROXIMATION LOGIC ---
            % We test if the phase shift is a rational fraction of Pi.
            
            % 1. Try theoretical N-fraction first (k/N * pi) - The Strong Hypothesis
            k_N = round(raw_phase * N / pi);
            err_N = rad2deg(abs(raw_phase - k_N * pi / N));
            
            if err_N < 2.5
                % Good fit with N-scaling
                fit_error_deg = err_N;
                if k_N == 0, n_frac_str = '0';
                elseif k_N == N, n_frac_str = 'pi';
                else, n_frac_str = sprintf('%d/%d pi', k_N, N);
                end
            else
                % 2. Try best rational fit (p/q * pi) using continued fractions
                % Target tolerance: 0.5 degrees relative to 180 (pi)
                tol = 2.5 / 180; 
                [num, den] = rat(raw_phase / pi, tol);
                
                err_rat = rad2deg(abs(raw_phase - (num/den)*pi));
                
                if err_rat < 2.5
                    fit_error_deg = err_rat;
                    if num == 0, n_frac_str = '0';
                    elseif num == den, n_frac_str = 'pi';
                    elseif den == 1, n_frac_str = sprintf('%d pi', num);
                    else, n_frac_str = sprintf('%d/%d pi', num, den);
                    end
                else
                    % No simple fraction found within tolerance
                    n_frac_str = '-';
                    fit_error_deg = err_N; % Report deviation from N-law for context
                end
            end
            % ---------------------------------------
            
            % Print formatted row
            fprintf('|  D%02d  |  %+2d   |   %6.1f    |   %5.2f   |   %-11s |     %5.3f     |\n', ...
                k, slope_val, deg_val, iso_spread_deg, n_frac_str, fit_error_deg);
            
            % Store for saving
            IsomerParameters(end+1).BranchID = b; %#ok<AGROW>
            IsomerParameters(end).RotorID = k;
            IsomerParameters(end).Chirality = slope_val;
            IsomerParameters(end).PhaseShift_Rad = raw_phase;
            IsomerParameters(end).PhaseShift_Deg = deg_val;
            IsomerParameters(end).ClusterRMS_Deg = iso_spread_deg;
            IsomerParameters(end).N_Fraction_Str = n_frac_str;
            IsomerParameters(end).Fit_Error_Deg = fit_error_deg;
            
            % Add to table arrays
            Tab_Branch = [Tab_Branch; b]; %#ok<AGROW>
            Tab_Rotor = [Tab_Rotor; k]; %#ok<AGROW>
            Tab_Slope = [Tab_Slope; slope_val]; %#ok<AGROW>
            Tab_Deg = [Tab_Deg; deg_val]; %#ok<AGROW>
            Tab_Spread = [Tab_Spread; iso_spread_deg]; %#ok<AGROW>
            Tab_NFrac{end+1,1} = n_frac_str; %#ok<AGROW>
            Tab_FitErr = [Tab_FitErr; fit_error_deg]; %#ok<AGROW>
        end
    end
    
    % Create Matlab Table object
    T_Isomers = table(Tab_Branch, Tab_Rotor, Tab_Slope, Tab_Deg, Tab_Spread, Tab_NFrac, Tab_FitErr, ...
        'VariableNames', {'BranchID', 'RotorID', 'Chirality', 'Phase_Deg', 'IsoSpread_Deg', 'Approx_N_Frac', 'Fit_Error_Deg'});
    
    % Save to .mat file
    save(dataFile, 'IsomerParameters', 'T_Isomers', '-append');
    fprintf('\n[INFO] Parameters saved to %s.\n', dataFile);

    %% 5. Visualization
    branchColors = lines(max(bestK, 4));
    
    if use3D
        plot3DTorus(N, selectedMasters, AngularCoords, Slopes, ClusterCentroids, bestK, branchColors);
    else
        plot2DSquare(N, selectedMasters, AngularCoords, Slopes, ClusterCentroids, bestK, branchColors);
    end
end

%% ========================================================================
%  ANALYSIS SUB-ROUTINES
%  ========================================================================

function Slopes = determineSlopes(AngularCoords, masterRef, N)
    % Uses Principal Component Analysis (PCA) on local neighborhoods of the 
    % manifold data to determine the tangent vector (Slope/Chirality).
    % Slope is either +1 (Co-rotating) or -1 (Counter-rotating).
    
    fprintf('Probing manifold slopes...\n');
    FullCircleCoords = 2 * AngularCoords; 
    M_total = size(FullCircleCoords, 1);
    
    n_probes = 50; 
    k_neighbors = 30;
    slope_votes = zeros(n_probes, N);
    
    for p = 1:n_probes
        % Pick a random point and find its neighbors
        idx_center = randi(M_total);
        center_pt = FullCircleCoords(idx_center, :);
        
        dists = zeros(M_total, 1);
        for j = 1:M_total
            diffs = abs(FullCircleCoords(j,:) - center_pt);
            diffs = min(diffs, 2*pi - diffs);
            dists(j) = sum(diffs.^2);
        end
        [~, sorted_idx] = sort(dists);
        neighbor_indices = sorted_idx(1:k_neighbors);
        
        % Perform PCA on local cloud
        local_cloud = FullCircleCoords(neighbor_indices, :);
        diff_cloud = mod(local_cloud - center_pt + pi, 2*pi) - pi;
        [coeff, ~, ~] = pca(diff_cloud);
        tangent = coeff(:, 1);
        
        % Normalize direction relative to Master
        if tangent(masterRef) < 0, tangent = -tangent; end
        slope_votes(p, :) = sign(tangent)';
    end
    % Take the mode (most common vote) as the true slope
    Slopes = mode(slope_votes, 1);
end

function [bestK, ClusterIdx, ClusterCentroids] = performClustering(RelativePhases, masterRef, n_kept, N)
    % Performs K-Means clustering in the Relative Phase Space to find Isomers.
    % Uses Silhouette Scores to determine the optimal number of clusters (K).
    
    eval_cols = setdiff(1:N, masterRef);
    DataForClustering = RelativePhases(:, eval_cols);
    % Map to Sine/Cosine to handle circular periodicity during clustering
    DataMetric = [cos(2*DataForClustering), sin(2*DataForClustering)];
    
    maxK = min(8, n_kept - 1); 
    sil_scores = zeros(1, maxK);
    
    % Evaluate K = 2 to MaxK
    for k = 2:maxK
        try
            idx = kmeans(DataMetric, k, 'Replicates', 5, 'Display', 'off');
            s = silhouette(DataMetric, idx);
            sil_scores(k) = mean(s);
        catch, sil_scores(k) = -1; end
    end
    [bestScore, bestK] = max(sil_scores);
    
    % Post-processing: Check if clusters are too close (Artifacts)
    if bestK > 1
        [ClusterIdx, C_metric] = kmeans(DataMetric, bestK, 'Replicates', 5);
        MIN_SEP = 15.0; 
        min_dist = 2 * sin(deg2rad(MIN_SEP));
        if min(pdist(C_metric)) < min_dist
            fprintf('  -> Merging clusters (too close).\n');
            bestK = 1; ClusterIdx = ones(n_kept, 1);
        end
    else
        ClusterIdx = ones(n_kept, 1);
    end
    
    % Fallback if structure is too weak
    if bestScore < 0.35 && bestK > 1
        fprintf('  -> Weak structure (Score %.2f). Reverting to K=1.\n', bestScore);
        bestK = 1; ClusterIdx = ones(n_kept, 1);
    end
    
    % Compute Centroids (Circular Mean)
    ClusterCentroids = zeros(bestK, N);
    for k = 1:bestK
        mask = (ClusterIdx == k);
        points = RelativePhases(mask, :);
        z = exp(1i * 2 * points);
        mean_z = mean(z, 1);
        ClusterCentroids(k, :) = mod(angle(mean_z) / 2, pi);
    end
end

%% ========================================================================
%  PLOTTING ROUTINES
%  ========================================================================

function plot2DSquare(N, masters, AngularCoords, Slopes, Centroids, K, colors)
    % Visualization: Flat Torus (Square) representation
    % Shows the linear correlation lines wrapping around the edges.
    figure('Color', 'w', 'Name', '2D Phase Locking Correlations');
    tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');
    
    for mID = masters
        for oID = 1:N
            nexttile; hold on;
            % Plot Raw Data Points
            scatter(AngularCoords(:, oID), AngularCoords(:, mID), ...
                'SizeData', 20, ...
                'MarkerEdgeColor', [0.5 0.5 0.5], ...
                'MarkerFaceColor', [0.5 0.5 0.5], ...
                'MarkerFaceAlpha', 0.4);
            
            % Plot Theoretical Isomer Lines
            slope_rel = Slopes(oID) / Slopes(mID);
            for b = 1:K
                C_m = Centroids(b, mID); C_o = Centroids(b, oID);
                intercept = mod(C_o - slope_rel * C_m, pi);
                
                tm = linspace(0, pi, 300);
                to = mod(slope_rel * tm + intercept, pi);
                jumps = find(abs(diff(to)) > pi/2); to(jumps) = NaN; % Fix wrap-around lines
                plot(to, tm, '-', 'Color', colors(b,:), 'LineWidth', 2);
            end
            xlim([0 pi]); ylim([0 pi]); axis square; grid on;
            title(sprintf('\\theta_{%d} vs \\theta_{%d}', mID, oID), ...
                'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'tex');
        end
    end
end

function plot3DTorus(N, masters, AngularCoords, Slopes, Centroids, K, colors)
    % Visualization: 3D Torus embedding
    % Wraps the data onto a donut shape to visually demonstrate the continuity
    % of the periodic boundaries.
    figure('Color', 'w', 'Name', '3D Toroidal Manifold');
    tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'none');
    
    R = 3; r = 1;
    [U, V] = meshgrid(linspace(0, 2*pi, 30), linspace(0, 2*pi, 15));
    Xb = (R + r*cos(V)) .* cos(U); 
    Yb = (R + r*cos(V)) .* sin(U); 
    Zb = r * sin(V);
    
    for mID = masters
        for oID = 1:N
            nexttile; hold on;
            % Draw Transparent Torus Surface
            surf(Xb, Yb, Zb, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.05, ...
                'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
            
            % Map Data Points to Torus Surface
            u = 2 * AngularCoords(:, mID); v = 2 * AngularCoords(:, oID);
            Xd = (R + r*cos(v)) .* cos(u); 
            Yd = (R + r*cos(v)) .* sin(u); 
            Zd = r * sin(v);
            scatter3(Xd, Yd, Zd, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.1);
            
            % Draw Isomer Lines
            slope_rel = Slopes(oID) / Slopes(mID);
            for b = 1:K
                C_m = Centroids(b, mID); C_o = Centroids(b, oID);
                intercept = mod(C_o - slope_rel * C_m, pi);
                tm = linspace(0, pi, 300);
                
                uv = 2 * tm; 
                vv = 2 * mod(slope_rel * tm + intercept, pi);
                
                Xc = (R + r*cos(vv)) .* cos(uv); 
                Yc = (R + r*cos(vv)) .* sin(uv); 
                Zc = r * sin(vv);
                plot3(Xc, Yc, Zc, '-', 'Color', colors(b,:), 'LineWidth', 3);
            end
            axis equal; axis off; view(-37.5, 30);
            title(sprintf('\\theta_{%d} vs \\theta_{%d}', mID, oID));
        end
    end
end