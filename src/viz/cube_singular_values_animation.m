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

function cube_singular_values_animation()
% =========================================================================
% OMNIDESIGN ANIMATION - CUBE SINGULAR VALUE ANALYZER
% =========================================================================
%
% PROJECT:      OmniDesign Optimizer Framework
% DESCRIPTION:
%   This interactive tool visualizes the performance of the Cube Chassis (N=8)
%   under different topological configurations (Isomers).
%
%   It creates a 3x6 grid layout to compare:
%   - Columns 1-3: Three distinct "Optimal" Isomers found by the optimizer.
%   - Column 4:    A Toroidal path with an incorrect phase offset (Suboptimal).
%   - Columns 5-6: Completely random geometric configurations.
%
%   For each column, it displays:
%   1. ROW 1 (Kinematics): The physical rotor orientations in real-time.
%   2. ROW 2 (Ellipsoids): The Force (Cyan) and Moment (Magenta) manipulability 
%      ellipsoids derived from the Jacobian.
%   3. ROW 3 (Spectrum): A bar chart of the 6 Singular Values (Sigma).
%
% METRICS COMPUTED:
%   - Kappa (Condition Number): Max(S)/Min(S). Lower is better (1.0 is isotropic).
%   - LogVol (Log-Volume): Sum of log singular values. Higher is better.
%
% CONTROLS:
%   - Slider: Manually scrub through the phase angle (Lambda).
%   - Play/Stop: Auto-rotate the rotors.
%   - Speed Slider: Adjust animation playback speed.
% =========================================================================

    clear; close all; clc;
    
    % --- CONFIGURATION --------------------------------------------------
    S.camRotFactor = 0.15;    % Speed of camera orbit relative to lambda
    S.trailHorizon = 2*pi;    % Length of the motion trail (radians)
    S.trailRes     = 20;      % Resolution of the trail lines
    
    % PARAMETER: Arrow Length (Visual scale of force vectors)
    S.ArrowLen     = 0.5;     
    
    % 1. VISUAL LAYOUT GEOMETRY (Proportions of the grid)
    RowRatios = [1.3, 1.3, 0.6]; 
    Gap_X = 0.0015;  
    Gap_Y = -0.03; 
    Row2_Lift = 0.03;      
    
    % Margins & Alignment
    Marg_Bot = 0.1; 
    Marg_Top = 0.04; 
    LabelColW = 0.04;  
    Marg_L   = 0.01; 
    Marg_R   = 0.02; 
    
    % 2. FONT SIZES
    FS.Title    = 14;  
    FS.RowLbl   = 16;  
    FS.Metrics  = 15;   
    FS.Axis     = 16;
    FS.momforce = 14;
    
    % 3. LOGIC THRESHOLDS
    MetricTol = 0.05; % Tolerance to highlight "Best" metric in green
    
    % --- Simulation State ---
    S.lambda = 0;
    S.isRunning = false;
    S.tLast = 0;
    
    % --- Geometry Setup (Cube) ---
    % Vertices of a cube inscribed in a sphere of radius 1 (normalized)
    S.P = [-1, 1, -1, 1, -1, 1, -1, 1; 1, 1, 1, 1, -1, -1, -1, -1; 1, 1, -1, -1, 1, 1, -1, -1] / sqrt(3); 
    N = 8;
    
    % --- Define Local Tangent Basis (U, V) ---
    % Defines the tangent plane at each vertex where the rotor force vector creates
    % a circle. U and V are orthogonal unit vectors tangent to the sphere.
    S.Basis_U = zeros(3, N); S.Basis_V = zeros(3, N);
    for k = 1:N
        n = S.P(:, k) / norm(S.P(:, k)); 
        if abs(dot(n, [0;0;1])) > 0.9, temp_ref = [1;0;0]; else, temp_ref = [0;0;1]; end
        u = cross(temp_ref, n); u = u / norm(u); v = cross(n, u);
        S.Basis_U(:, k) = u; S.Basis_V(:, k) = v;
    end
    
    % --- Phase Offsets & Random Generation ---
    % Hardcoded phase offsets for the specific Isomers (Columns 1-4)
    S.Offsets = zeros(N, 6);
    S.Offsets(:, 1) = [0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5]' * pi;           % Isomer 1
    S.Offsets(:, 2) = [0, 2/3, 2/3, 0, 1/6, 0.5, 0.5, 1/6]' * pi;       % Isomer 2
    S.Offsets(:, 3) = [0, 5/6, 1/3, 0.5, 1/3, 0.5, 0, 5/6]' * pi;       % Isomer 3
    vals = linspace(0, 2*pi, N+1)'; S.Offsets(:, 4) = vals(1:end-1);    % Torus Path (Bad Offset)
    
    rng(101); 
    S.RandPhase1 = rand(N,1) * 2*pi; % Phasing for Col 5
    S.RandPhase2 = rand(N,1) * 2*pi; % Phasing for Col 6
    
    % --- NEW: RANDOM 3D BASIS FOR COLUMN 5 (Arbitrary 3D Rotation) ---
    % Unlike cols 1-4 which rotate tangent to the chassis surface, 
    % these columns simulate completely broken/randomized rotor mounts.
    S.Rand3D_U = zeros(3, N);
    S.Rand3D_V = zeros(3, N);
    for k = 1:N
        % Random normal vector (axis of rotation)
        axisRot = randn(3,1); axisRot = axisRot/norm(axisRot);
        % Create two orthogonal vectors in the plane perpendicular to axisRot
        if abs(dot(axisRot, [0;0;1])) > 0.9, temp = [1;0;0]; else, temp = [0;0;1]; end
        u_r = cross(axisRot, temp); u_r = u_r/norm(u_r);
        v_r = cross(axisRot, u_r);
        S.Rand3D_U(:, k) = u_r;
        S.Rand3D_V(:, k) = v_r;
    end
    % --- NEW: RANDOM 3D BASIS FOR COLUMN 6 (Arbitrary 3D Rotation) ---
    S.Rand3D_U2 = zeros(3, N);
    S.Rand3D_V2 = zeros(3, N);
    for k = 1:N
        % Random normal vector (axis of rotation)
        axisRot = randn(3,1); axisRot = axisRot/norm(axisRot);
        % Create two orthogonal vectors in the plane perpendicular to axisRot
        if abs(dot(axisRot, [0;0;1])) > 0.9, temp = [1;0;0]; else, temp = [0;0;1]; end
        u_r = cross(axisRot, temp); u_r = u_r/norm(u_r);
        v_r = cross(axisRot, u_r);
        S.Rand3D_U2(:, k) = u_r;
        S.Rand3D_V2(:, k) = v_r;
    end
    
    % --- UI Setup ---
    fig = figure('Name', 'Cube Analyzer v15', 'Color', 'w', ...
                 'Position', [10, 50, 1800, 900], 'NumberTitle', 'off');
             
    colTitles = {'1. OPTIMAL ISOMER 1', '2. OPTIMAL ISOMER 2', '3. OPTIMAL ISOMER 3', ...
                 '4. ON TORUS w/ WRONG OFFSET ', '5. RANDOM 1', '6. RANDOM 2'};
    
    % Object Handles
    H.axChassis = gobjects(1,6); H.axEllips = gobjects(1,6); H.axHist = gobjects(1,6);
    G.arrows = gobjects(6, N); G.lines = gobjects(6, N); G.trails = gobjects(6, N);
    G.ellipForce = gobjects(6, 1); G.ellipMoment = gobjects(6, 1);
    G.bars = gobjects(6, 1); G.txtMetrics = gobjects(6, 1);
    [Ex, Ey, Ez] = sphere(20); % Base sphere for ellipsoids
    
    % --- Calculate Grid Geometry ---
    TotalAvailableH = 1 - Marg_Bot - Marg_Top - (2 * Gap_Y);
    NormRowH = (RowRatios / sum(RowRatios)) * TotalAvailableH;
    ColW = (1 - Marg_L - LabelColW - Marg_R - (5 * Gap_X)) / 6;
    
    % --- Create Label Column (Row Labels) ---
    commonX = Marg_L + (LabelColW/2);
    annotation('textarrow',[commonX commonX],[0.85 0.85], 'String','KINEMATICS', ...
        'HeadStyle','none','LineStyle','none','FontSize',FS.RowLbl,'FontWeight','bold','TextRotation',90, ...
        'HorizontalAlignment','center');
    annotation('textarrow',[commonX commonX],[0.50 0.50], 'String','3D PROJ. ELLIPSOIDS', ...
        'HeadStyle','none','LineStyle','none','FontSize',FS.RowLbl,'FontWeight','bold','TextRotation',90, ...
        'HorizontalAlignment','center');
    annotation('textarrow',[commonX commonX],[0.18 0.18], 'String','SINGULAR VALUES', ...
        'HeadStyle','none','LineStyle','none','FontSize',FS.RowLbl,'FontWeight','bold','TextRotation',90, ...
        'HorizontalAlignment','center');
        
    % --- Create Plots Loop ---
    for col = 1:6
        pos_x = Marg_L + LabelColW + (col-1) * (ColW + Gap_X);
        
        % ROW 1: KINEMATICS (3D Robot View)
        pos_y_1 = Marg_Bot + NormRowH(3) + Gap_Y + NormRowH(2) + Gap_Y;
        H.axChassis(col) = axes('Position', [pos_x, pos_y_1, ColW, NormRowH(1)]);
        hold on; 
        
        % Camera and Axis Properties
        axis equal; 
        set(gca, 'DataAspectRatio', [1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
        xlim([-2.0 2.0]); ylim([-2.0 2.0]); zlim([-2.0 2.0]); 
        view(45, 30); 
        camva(5); % Fixed FOV (prevent pulsating)
        axis off; 
        
        % Title Placement
        text(0, 0, 2.0, colTitles{col}, ...
            'HorizontalAlignment', 'center', ...
            'FontSize', FS.Title, ...
            'FontWeight', 'bold');
            
        % Draw Convex Hull (Chassis Body)
        k_hull = convhull(S.P(1,:), S.P(2,:), S.P(3,:));
        trisurf(k_hull, S.P(1,:), S.P(2,:), S.P(3,:), 'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'k');
        
        % Initialize Rotors (Arrows and Trails)
        for k=1:N
            % STYLING: Propeller 1 Highlight to track rotation
            if k == 1
                col_arrow = [0.6 0 0.8]; % Purple
                lw_arrow = 2.5;
                lw_trail = 2;
                col_trail = [0.6 0 0.8];
            else
                col_arrow = 'r';
                lw_arrow = 2.0;
                lw_trail = 2;
                col_trail = [0 0.4 0.9];
            end
            
            G.trails(col,k) = plot3(0,0,0, ':', 'Color', col_trail, 'LineWidth', lw_trail); 
            G.arrows(col,k) = quiver3(0,0,0, 0,0,0, 0, 'Color', col_arrow, 'LineWidth', lw_arrow, 'MaxHeadSize', 0.6, 'Clipping', 'off', 'AutoScale', 'off');
            G.lines(col,k) = plot3([0 0], [0 0], [0 0], '--', 'Color', col_arrow, 'LineWidth', 1.0);
        end
        
        % ROW 2: ELLIPSOIDS (Force and Moment Capability)
        pos_y_2 = Marg_Bot + NormRowH(3) + Gap_Y + Row2_Lift; 
        H.axEllips(col) = axes('Position', [pos_x, pos_y_2, ColW, NormRowH(2)]);
        hold on; axis equal; view(0, 0); xlim([-2.2 2.2]); ylim([-2.2 2.2]); zlim([-2.5 2]); axis off;
        
        G.ellipForce(col) = surf(Ex, Ey, Ez, 'FaceColor', 'c', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        G.ellipMoment(col) = surf(Ex, Ey, Ez, 'FaceColor', 'm', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        
        % Labels
        text(1.2, 0, 1.4, 'Force', 'FontWeight', 'bold', 'Color', 'k','FontSize', FS.momforce );
        text(1.2, 0, -1.4, 'Moment', 'FontWeight', 'bold', 'Color', 'k','FontSize', FS.momforce );
        
        lighting gouraud; camlight;
        
        % ROW 3: SPECTRUM (SVD Bar Chart)
        pos_y_3 = Marg_Bot;
        SlimW = ColW * 0.7; SlimX = pos_x + (ColW - SlimW) / 2; 
        H.axHist(col) = axes('Position', [SlimX, pos_y_3, SlimW, NormRowH(3)]);
        hold on; grid on; ylim([0 4.0]); xlim([0.5 6.5]); set(gca, 'XTick', 1:6, 'FontSize', FS.Axis);
        
        % Restore Y-Label on first plot only
        if col == 1
            ylabel('\sigma Value', 'FontSize', FS.Axis);
        else
            set(gca, 'YTickLabel', []); 
        end
        
        G.bars(col) = bar(1:6, zeros(1,6), 'FaceColor', [0.2 0.6 0.8], 'BarWidth', 0.5);
        G.txtMetrics(col) = text(1.0, 3.5, '', 'FontSize', FS.Metrics, 'BackgroundColor', 'w');
    end
    
    % --- UI CONTROLS ---
    hSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 2*pi, 'Value', 0, ...
        'Position', [600, 25, 400, 20], 'Callback', @slider_cb); 
    
    hSpeed = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 0.15, 'Value', 0.04, ...
        'Position', [1050, 25, 100, 20]);
    
    uicontrol('Style', 'pushbutton', 'String', 'PLAY', 'Position', [100, 20, 80, 40], ...
        'FontSize', FS.RowLbl, 'Callback', @play_cb);
    uicontrol('Style', 'pushbutton', 'String', 'STOP', 'Position', [200, 20, 80, 40], ...
        'FontSize', FS.RowLbl, 'Callback', @pause_cb);
        
    % --- CORE UPDATE LOGIC ---
    function update_viz(lambda)
        S.lambda = lambda;
        shift_dist_z = 1; 
        
        % Camera Orbit Logic (Top Row Only) - slowly rotates view
        camAngle = 45 + rad2deg(lambda * S.camRotFactor);
        
        AllLogVol = zeros(1,6);
        AllKappa  = zeros(1,6); 
        
        for c = 1:6
            set(H.axChassis(c), 'View', [camAngle, 30]);
            
            CurrentBasis = zeros(3, N);
            
            for k = 1:N
                % 1. DEFINE BASIS (u, v) AND OFFSET FOR THIS COLUMN
                u_base = S.Basis_U(:,k); 
                v_base = S.Basis_V(:,k);
                off = 0;
                
                if c <= 4
                    % Standard Tangent Rotation (Isomers)
                    off = S.Offsets(k, c);
                    u = u_base; v = v_base;
                elseif c == 5
                    % RANDOM 3D ROTATION 1 (Forces rotate about random 3D axes)
                    off = S.RandPhase1(k);
                    u = S.Rand3D_U(:,k); 
                    v = S.Rand3D_V(:,k); 
                elseif c == 6
                    % RANDOM 3D ROTATION 2 (Forces rotate about random 3D axes)
                    off = S.RandPhase2(k);
                    u = S.Rand3D_U2(:,k); 
                    v = S.Rand3D_V2(:,k); 
                end
                
                % 2. TRAIL CALCULATION (History)
                % Draws the "tail" behind the vector to visualize the manifold curve
                t_hist = linspace(lambda - S.trailHorizon, lambda, S.trailRes);
                trailPts = zeros(3, S.trailRes);
                for i = 1:S.trailRes
                    theta_h = t_hist(i) + off;
                    vec_h = cos(theta_h)*u + sin(theta_h)*v;
                    trailPts(:,i) = S.P(:,k) + vec_h * S.ArrowLen;
                end
                set(G.trails(c,k), 'XData', trailPts(1,:), 'YData', trailPts(2,:), 'ZData', trailPts(3,:));
                
                % 3. CURRENT VECTOR UPDATE
                theta = lambda + off;
                vec = cos(theta)*u + sin(theta)*v;
                CurrentBasis(:,k) = vec;
                
                % Update Quiver Arrows
                set(G.arrows(c,k), 'XData', S.P(1,k), 'YData', S.P(2,k), 'ZData', S.P(3,k), ...
                    'UData', vec(1)*S.ArrowLen, 'VData', vec(2)*S.ArrowLen, 'WData', vec(3)*S.ArrowLen);
                
                set(G.lines(c,k), 'XData', [S.P(1,k) S.P(1,k)-0.5*vec(1)*S.ArrowLen], ...
                                  'YData', [S.P(2,k) S.P(2,k)-0.5*vec(2)*S.ArrowLen], ...
                                  'ZData', [S.P(3,k) S.P(3,k)-0.5*vec(3)*S.ArrowLen]);
            end
            
            % 4. BUILD MATRIX A (The Grasp Matrix)
            % Maps rotor forces (u) to Body Wrench (Force + Moment)
            A = zeros(6, N);
            for k=1:N
                f = CurrentBasis(:,k); 
                m = cross(S.P(:,k), f);
                A(1:3, k) = f; 
                A(4:6, k) = m;
            end
            
            % 5. SVD & METRICS
            sv = svd(A); 
            AllLogVol(c) = sum(log(sv + 1e-9)); % Log-Volume (Proposed Metric)
            AllKappa(c)  = max(sv)/min(sv);     % Condition Number (Standard Metric)
            
            % Update Ellipsoids (Visualization of A*A')
            % Force Ellipsoid
            G_f = A(1:3,:) * A(1:3,:)'; [V, D] = eig(G_f); r = sqrt(diag(D)) * 0.45; 
            pts = V * diag(r) * [Ex(:)'; Ey(:)'; Ez(:)'];
            set(G.ellipForce(c), 'XData', reshape(pts(1,:), size(Ex)), ...
                                 'YData', reshape(pts(2,:), size(Ey)), ...
                                 'ZData', reshape(pts(3,:), size(Ez)) + shift_dist_z); 
            
            % Moment Ellipsoid
            G_m = A(4:6,:) * A(4:6,:)'; [V, D] = eig(G_m); r = sqrt(diag(D)) * 0.45;
            pts = V * diag(r) * [Ex(:)'; Ey(:)'; Ez(:)'];
            set(G.ellipMoment(c), 'XData', reshape(pts(1,:), size(Ex)), ...
                                  'YData', reshape(pts(2,:), size(Ey)), ...
                                  'ZData', reshape(pts(3,:), size(Ez)) - shift_dist_z); 
            
            % Update Histogram Bars
            set(G.bars(c), 'YData', sv);
        end
        
        % Highlight the "Best" result dynamically
        BestVol = max(AllLogVol);
        for c = 1:6
            isBest = (AllLogVol(c) >= BestVol - MetricTol);
            if isBest, colTxt = [0 0.6 0]; w = 'bold'; else, colTxt = [0.8 0 0]; w='normal'; end
            set(G.txtMetrics(c), 'String', sprintf('\\kappa=%.1f\nLogVol=%.1f', AllKappa(c), AllLogVol(c)), ...
                'Color', colTxt, 'FontWeight', w);
        end
        drawnow limitrate;
    end
    
    function pause_cb(~, ~), S.isRunning = false; end
    
    % Slider Callback (Manual Scrubbing)
    function slider_cb(src, ~)
        S.isRunning = false; 
        val = get(src, 'Value');
        update_viz(val);
    end
    
    % Play Callback (Animation Loop)
    function play_cb(~, ~)
        if S.isRunning, return; end
        S.isRunning = true; S.tLast = tic;
        while S.isRunning
            if ~ishandle(fig), S.isRunning = false; break; end
            freq = get(hSpeed, 'Value');
            dt = toc(S.tLast); S.tLast = tic;
            val = mod(get(hSlider, 'Value') + (freq * pi * dt), 2*pi);
            set(hSlider, 'Value', val);
            update_viz(val);
            pause(0.01);
        end
    end
    
    % Initialize
    update_viz(0);
end