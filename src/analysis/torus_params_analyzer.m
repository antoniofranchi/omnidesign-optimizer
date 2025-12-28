function OmniDesign_Analysis_TorusParams_intrinsic()
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
% OMNIDESIGN OPTIMIZER - MODULE 2: MANIFOLD PARAMETERIZATION
% =========================================================================
%
% PROJECT:      OmniDesign Optimizer Framework
% DESCRIPTION:
%   This script bridges the gap between the raw optimization results (Module 1)
%   and the topological analysis (Module 3). It transitions the data from 
%   the extrinsic 2D unit disc representation to the intrinsic angular 
%   coordinates of the Tangent Torus (T^N).
%
%   Unlike statistical approaches, this module constructs the EXACT 
%   geometric tangent basis {u_i, v_i} for every rotor location P_i defined 
%   by the chassis geometry.
%
% METHODOLOGY:
%   1. LOADING: Imports the pruned 'PointsTable' and chassis 'P'.
%   2. BASIS CONSTRUCTION: Defines a local orthonormal frame for each rotor:
%      - n: Normal vector (Radial direction from center).
%      - u: Tangent vector parallel to the global XY plane (Azimuthal).
%      - v: Tangent vector orthogonal to u and n (Poloidal).
%   3. PROJECTION: Projects every optimal force vector d_i onto this basis.
%   4. PARAMETERIZATION: Computes the intrinsic phase angle:
%      theta_i = atan2( d . v, d . u )
%
% INPUT:
%   - 'PointsTable' and 'P' (from Module 1 .mat file).
%
% OUTPUT:
%   - 'TorusMapping_Data.mat': Contains 'AngularCoords' (the intrinsic theta 
%     values) required for the Topological Analysis.
%
% =========================================================================

    %% 1. Data Loading
    %  Checks workspace for existing data; otherwise prompts user to load.
    if evalin('base', 'exist(''PointsTable'',''var'')') && evalin('base', 'exist(''P'',''var'')')
        PointsTable = evalin('base', 'PointsTable');
        P = evalin('base', 'P');
        disp('Loaded data from workspace.');
    else
        [file, path] = uigetfile('*.mat', 'Select Optimizer Results File');
        if isequal(file,0), disp('Selection cancelled.'); return; end
        loadedData = load(fullfile(path, file));
        
        if isfield(loadedData, 'PointsTable') && isfield(loadedData, 'P')
            PointsTable = loadedData.PointsTable;
            P = loadedData.P;
        else
            error('Invalid file: Must contain ''PointsTable'' and ''P''.');
        end
        disp(['Loaded ' file]);
    end
    
    N = size(P, 2);
    fprintf('Detected N = %d propellers.\n', N);

    %% 2. User Configuration (Data Pruning)
    %  Interactive dialog to filter the dataset. This allows the user to
    %  focus only on the highest-quality solutions (the "Global Manifold")
    %  by removing suboptimal points based on Reward (Volume) and Cost.
    d = dialog('Position',[300 300 300 150], 'Name', 'Step 1: Data Pruning');
    uicontrol('Parent',d,'Style','text','Position',[20 100 250 20],...
        'String','Keep Top % (Reward / Cost):','HorizontalAlignment','left');
    edit_reward = uicontrol('Parent',d,'Style','edit','Position',[20 80 50 20],'String','100');
    edit_cost   = uicontrol('Parent',d,'Style','edit','Position',[80 80 50 20],'String','100');
    uicontrol('Parent',d,'Position',[80 20 140 30],'String','Run Projection',...
        'Callback','uiresume(gcbf)');
    uiwait(d);
    
    if ~ishandle(d), disp('Cancelled.'); return; end
    pct_reward = str2double(get(edit_reward,'String'));
    pct_cost   = str2double(get(edit_cost,'String'));
    delete(d);

    %% 3. Filter Data
    %  Apply the percentiles defined above to mask the data.
    RawMetrics = table2array(PointsTable(:, end-1:end)); 
    thresh_reward = prctile(RawMetrics(:, 1), 100 - pct_reward);
    thresh_cost   = prctile(RawMetrics(:, 2), pct_cost);
    
    valid_mask = (RawMetrics(:, 1) >= thresh_reward) & (RawMetrics(:, 2) <= thresh_cost);
    
    % Extract only the coordinate columns for processing
    varNames = PointsTable.Properties.VariableNames;
    isCoordCol = contains(varNames, 'Prop') & (contains(varNames, '_x') | contains(varNames, '_y'));
    FilteredMatrix = table2array(PointsTable(valid_mask, isCoordCol));
    
    n_kept = size(FilteredMatrix, 1);
    if n_kept == 0, errordlg('No data points left after filtering!'); return; end
    fprintf('Filtering: Kept %d / %d solutions.\n', n_kept, height(PointsTable));

    %% 4. Geometric Basis Construction & Projection
    %  This loop computes the intrinsic angle Theta_i for every rotor i.
    %  It constructs the Tangent Space T_p S^2 exactly based on geometry.
    
    AngularCoords = zeros(n_kept, N);
    
    % Prepare storage structure for visualization frames
    BasisFrames = struct('u', cell(1,N), 'v', cell(1,N), 'n', cell(1,N));
    
    fprintf('Projecting data onto Tangent Torus frames...\n');
    
    for i = 1:N
        % --- A. Construct Local Basis {u, v} ---
        
        % 1. Normal vector (n): The radial direction of the vertex
        n_vec = P(:, i) / norm(P(:, i));
        
        % 2. Tangent vector (u): Defined parallel to Global XY plane.
        %    Calculated via cross product with Global Z (k_hat).
        k_hat = [0; 0; 1];
        
        % Handle Singularity: If point is exactly at the pole (0,0,1)
        if abs(dot(n_vec, k_hat)) > 0.99
            u_vec = [1; 0; 0]; % Arbitrary valid tangent for poles
        else
            u_vec = cross(k_hat, n_vec);
            u_vec = u_vec / norm(u_vec);
        end
        
        % 3. Tangent vector (v): Completes the orthonormal set ("Tangent Up")
        v_vec = cross(n_vec, u_vec);
        v_vec = v_vec / norm(v_vec);
        
        % Store basis for later 3D visualization
        BasisFrames(i).u = u_vec;
        BasisFrames(i).v = v_vec;
        BasisFrames(i).n = n_vec;
        
        % --- B. Reconstruct & Project Data ---
        col_x = 2*i - 1;
        col_y = 2*i;
        
        % Retrieve Extrinsic 2D disk coordinates
        x_disk = FilteredMatrix(:, col_x);
        y_disk = FilteredMatrix(:, col_y);
        
        % Reconstruct full 3D vector d (Canonical Mapping z > 0)
        % z = sqrt(1 - x^2 - y^2)
        r2 = x_disk.^2 + y_disk.^2;
        r2(r2 > 1) = 1; % Numerical safety for precision errors
        z_disk = sqrt(1 - r2);
        
        D_matrix = [x_disk, y_disk, z_disk]'; % 3 x M
        
        % Project onto local basis components (Dot Products)
        % coord_u = d . u
        % coord_v = d . v
        comp_u = u_vec' * D_matrix; % 1 x M
        comp_v = v_vec' * D_matrix; % 1 x M
        
        % Compute Intrinsic Phase Angle
        theta = atan2(comp_v, comp_u); % Result in (-pi, pi]
        
        % --- C. Map to RP1 (Projective Line) ---
        % Since line-of-force d is equivalent to -d, the domain is [0, pi).
        % We standardize the angle to the upper semi-circle.
        theta(theta < 0) = theta(theta < 0) + pi;
        
        AngularCoords(:, i) = theta(:);
    end
    
    %% 5. Visualization 1: Intrinsic 2D Projection
    %  Plots the distribution of angles on the abstract S1 circle.
    plotTangentProjection(N, AngularCoords, BasisFrames);
    
    %% 6. Visualization 2: 3D Geometry & Clouds
    %  Plots the 3D chassis with the tangent frames and solution clouds.
    plot3DGeometryAndClouds(P, BasisFrames, AngularCoords, N);

    %% 7. Save Data
    saveName = 'TorusMapping_Data.mat';
    % We save P and BasisFrames so Module 3 knows the topology
    save(saveName, 'AngularCoords', 'BasisFrames', 'FilteredMatrix', 'P', 'N', 'n_kept');
    
    fprintf('\nSUCCESS: Analysis Part 1 Complete.\n');
    fprintf(' - Intrinsic coordinates extracted using exact geometric frames.\n');
    fprintf(' - Data saved to ''%s''.\n', saveName);
    fprintf(' - Next Step: Run Module 3 (OmniDesign_Analysis_N5Law)\n');
end

%% ========================================================================
%  HELPER FUNCTIONS: VISUALIZATION
%  ========================================================================

function plotTangentProjection(N, AngularCoords, BasisFrames)
    % Visualize the distribution of angles on the intrinsic circle
    f = figure('Name', 'Tangent Torus Projection (2D)', 'Color', 'w', 'Position', [100 100 1200 600]);
    t = tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'loose');
    
    for i = 1:N
        nexttile; hold on; axis equal;
        
        % 1. Plot Unit Circle (The Tangent Space RP1)
        theta_circ = linspace(0, 2*pi, 200);
        plot(cos(theta_circ), sin(theta_circ), 'k-', 'Color', [0.8 0.8 0.8]);
        
        % 2. Plot Histogram of Angles on the Circle
        th = AngularCoords(:, i);
        
        % We plot both theta and theta+pi to show the full symmetric distribution
        % visually (handling the RP1 equivalence).
        u_pts = [cos(th); cos(th + pi)];
        v_pts = [sin(th); sin(th + pi)];
        
        scatter(u_pts, v_pts, 10, [0 0.44 0.74], 'filled', 'MarkerFaceAlpha', 0.1);
        
        % 3. Draw Basis Vectors
        quiver(0,0, 1,0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5); % u
        quiver(0,0, 0,1, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5); % v
        
        xlim([-1.2 1.2]); ylim([-1.2 1.2]); axis off;
        title(sprintf('Rotor %d (\\theta_i)', i), 'FontSize', 12);
    end
    sgtitle('Intrinsic Phase Distribution on Tangent Bundle', 'FontSize', 16, 'FontWeight', 'bold');
end

function plot3DGeometryAndClouds(P, BasisFrames, AngularCoords, N)
    % 3D Visualization of the chassis, local frames, and solution clouds
    fig = figure('Name', '3D Tangent Bundle & Chassis', 'Color', 'w', 'Position', [150 150 1000 800]);
    hold on; axis equal; grid on;
    view(45, 30);
    
    % --- 1. Draw Chassis Geometry ---
    % Center
    scatter3(0,0,0, 80, 'k', 'filled'); 
    
    for i = 1:N
        p_curr = P(:, i);
        p_next = P(:, mod(i, N) + 1);
        
        % Spoke (Center to Vertex)
        plot3([0, p_curr(1)], [0, p_curr(2)], [0, p_curr(3)], ...
              'k--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
          
        % Rim (Vertex to Vertex)
        plot3([p_curr(1), p_next(1)], [p_curr(2), p_next(2)], [p_curr(3), p_next(3)], ...
              'k-', 'LineWidth', 2);
    end
    
    % --- 2. Draw Frames and Clouds ---
    s_basis = 0.4;  % Length of basis vectors for viz
    s_cloud = 0.3;  % Radius of the data cloud ring
    
    for i = 1:N
        pt = P(:, i);
        u = BasisFrames(i).u;
        v = BasisFrames(i).v;
        n = BasisFrames(i).n;
        
        % A. Draw Basis Vectors
        % u (Red) - Tangent Horizontal
        quiver3(pt(1), pt(2), pt(3), u(1), u(2), u(3), s_basis, ...
            'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
            
        % v (Green) - Tangent Vertical
        quiver3(pt(1), pt(2), pt(3), v(1), v(2), v(3), s_basis, ...
            'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
            
        % n (Blue Dashed) - Discarded Radial/Normal
        quiver3(pt(1), pt(2), pt(3), n(1), n(2), n(3), s_basis, ...
            'b', 'LineWidth', 1, 'LineStyle', '--', 'MaxHeadSize', 0.3);
            
        % B. Labels
        % Label Vertex Theta
        text(pt(1)*1.2, pt(2)*1.2, pt(3)*1.2 + 0.1, ...
             sprintf('\\theta_{%d}', i), 'FontSize', 14, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k');
         
        % Label Vectors (Small offset)
        text(pt(1)+u(1)*s_basis, pt(2)+u(2)*s_basis, pt(3)+u(3)*s_basis, 'u', 'Color', 'r', 'FontSize', 8);
        text(pt(1)+v(1)*s_basis, pt(2)+v(2)*s_basis, pt(3)+v(3)*s_basis, 'v', 'Color', 'g', 'FontSize', 8);
        
        % C. Draw Data Clouds 
        % We project the theta values back into 3D using the local basis
        thetas = AngularCoords(:, i);
        
        % Downsample for rendering speed if necessary
        if length(thetas) > 2000
            idx = randperm(length(thetas), 2000);
            thetas = thetas(idx);
        end
        
        % Create circular distribution on the tangent plane
        % P_cloud = P_vertex + R * (cos(theta)*u + sin(theta)*v)
        
        % We plot theta AND theta+pi to show the full RP1 line symmetry
        % (First Half)
        cx1 = pt(1) + s_cloud * (cos(thetas).*u(1) + sin(thetas).*v(1));
        cy1 = pt(2) + s_cloud * (cos(thetas).*u(2) + sin(thetas).*v(2));
        cz1 = pt(3) + s_cloud * (cos(thetas).*u(3) + sin(thetas).*v(3));
        
        % (Second Half - Antipodal)
        cx2 = pt(1) + s_cloud * (cos(thetas+pi).*u(1) + sin(thetas+pi).*v(1));
        cy2 = pt(2) + s_cloud * (cos(thetas+pi).*u(2) + sin(thetas+pi).*v(2));
        cz2 = pt(3) + s_cloud * (cos(thetas+pi).*u(3) + sin(thetas+pi).*v(3));
        
        % Plot Points
        scatter3([cx1; cx2], [cy1; cy2], [cz1; cz2], 3, [0 0.44 0.74], ...
                 'filled', 'MarkerFaceAlpha', 0.3);
    end
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Geometric Context: Vertices, Tangent Frames & Solution Clouds');
    lighting gouraud;
end