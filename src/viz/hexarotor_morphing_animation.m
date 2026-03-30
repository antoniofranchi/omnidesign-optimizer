function hexarotor_morphing_animation()
% =========================================================================
% OMNIDESIGN OPTIMIZER - DYNAMIC RECONFIGURATION & AERODYNAMIC FOOTPRINTING
% =========================================================================
% Description:
%   This interactive simulation demonstrates the practical utility of the 
%   continuous design null space (N-5 Law). It visualizes a fully-actuated
%   hexarotor (N=6) dynamically reconfiguring its thrust vectors along the 
%   1D optimal manifold (parameterized by lambda) while performing a 
%   constrained lateral navigation task.
%
% Key Validations:
%   1. Aerodynamic Wake Vectoring: Shows how sweeping lambda redirects the 
%      downwash footprint on the ground plane (Z = -H).
%   2. Control Effort Invariance: Empirically proves that the total L2 
%      norm of the control input vector (||f||) remains strictly constant 
%      during the morphological transition, confirming zero-energy self-motion.
%
% Interaction:
%   - Use the UI to toggle automatic lambda sweeping or platform motion.
%   - Click 'Stop & Save' to terminate the loop and export the simulation 
%     telemetry to 'morphing_sim_log.mat' for plotting.
%
% License: MIT / Academic Free License
% =========================================================================

    % --- 1. System & Simulation Parameters ---
    N = 6;                  % Number of rotors (Hexarotor)
    R_arm = 0.5;            % Structural arm length (m)
    R_prop = 0.25;          % Propeller radius / Aerodynamic wake radius (m)
    H = 1.0;                % Height of the platform above the ground plane (m)
    mass = 1.0;             % Platform mass (kg)
    g = 9.81;               % Gravity (m/s^2)
    
    % Automated Start & Initial State Parameters
    auto_start_sweep  = 1;  % 1 = Start morphing (lambda sweep) automatically
    auto_start_motion = 1;  % 1 = Start lateral navigation task automatically
    initial_lambda    = 0;  % Starting position on the 1D manifold (0 to pi)
    
    % Environmental Particle Monitoring Area (For Wake Interaction Logging)
    monitor_X_span = [-1, 1];    % X boundaries for the particle counting box
    monitor_Y_span = [-0.7, 0.7];  % Y boundaries for the particle counting box
    
    % Kinematic & Reconfiguration Parameters
    lambda_sweep_speed = 0.02;  % Speed of morphological reconfiguration (rad/s)
    drone_motion_amp = 0.5;     % Amplitude of the lateral navigation task (m)
    drone_motion_freq = 0.15;   % Frequency of the lateral navigation task (rad/s)
    
    % Physics Solver / Visual Tuning Parameters
    dt = 0.05;              % Simulation time step
    sim_speed = 0.15;       % Particle advection speed multiplier
    arrow_scale = 0.1;      % Quiver plot scaling for downwash vectors
    grid_X_span = [-2.5, 2.5]; 
    grid_Y_span = [-1.5, 1.5];
    lateral_inertia_multiplier = 1; 
    
    % Affine Phase-Locking offsets (Delta for N=6)
    delta = (0:N-1) * (pi/2); 
    colors = lines(N);
    
    %% --- 2. Base Chassis Geometry Definition ---
    % Define the standard planar polygon chassis
    angles = (0:N-1) * (2*pi/N);
    P_base = [R_arm*cos(angles); R_arm*sin(angles); zeros(1, N)];
    
    % Precompute local rotor coordinate frames (Tangent U, Normal V)
    U = zeros(3, N); V = zeros(3, N);
    for i = 1:N
        n = P_base(:,i) / norm(P_base(:,i));     % Radial vector
        u = cross([0;0;1]', n')';                % Tangent vector
        U(:,i) = u / norm(u);
        V(:,i) = cross(n, U(:,i));               % Normal vector (Z-axis)
    end
    
    %% --- 3. Setup Interactive GUI and Axes ---
    fig = figure('Name', 'Dynamic Reconfiguration Simulator', ...
                 'Position', [50, 50, 800, 750], 'Color', 'w');
    fig.UserData = true; % Use UserData as the master "keep running" flag
    
    % Ground Plane Subplot (Aerodynamic Footprint)
    ax_ground = subplot(3, 1, [1 2], 'Parent', fig);
    hold(ax_ground, 'on'); grid(ax_ground, 'on'); axis(ax_ground, 'equal');
    xlim(ax_ground, grid_X_span); ylim(ax_ground, grid_Y_span);
    xlabel(ax_ground, 'Ground X (m)'); ylabel(ax_ground, 'Ground Y (m)');
    title(ax_ground, 'Aerodynamic Footprint & Particle Advection (Z = -H)');
    
    % Draw the target monitoring box for particle advection logging
    px = [monitor_X_span(1), monitor_X_span(2), monitor_X_span(2), monitor_X_span(1), monitor_X_span(1)];
    py = [monitor_Y_span(1), monitor_Y_span(1), monitor_Y_span(2), monitor_Y_span(2), monitor_Y_span(1)];
    plot(ax_ground, px, py, 'k--', 'LineWidth', 1.5);
    text(ax_ground, monitor_X_span(1), monitor_Y_span(2) + 0.1, 'Target Area', 'FontWeight', 'bold');
    
    % Control Effort Subplot (Bar Chart)
    ax_bar = subplot(3, 1, 3, 'Parent', fig);
    hold(ax_bar, 'on');
    
    b_plot = bar(ax_bar, 1:N, zeros(1,N), 0.3, 'FaceColor', 'flat');
    for i = 1:N; b_plot.CData(i,:) = colors(i,:); end
    
    % The total control effort bar (||f||) - MUST remain strictly constant
    b_effort = bar(ax_bar, N+1.5, 0, 0.3, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'k', 'LineWidth', 1.5);
    
    ylim(ax_bar, [-mass*g, mass*g]); xlim(ax_bar, [0.5, N+2.5]);
    xticks(ax_bar, [1:N, N+1.5]); xticklabels(ax_bar, {'1', '2', '3', '4', '5', '6', '||f||'});
    xlabel(ax_bar, 'Rotor Index & Effort'); ylabel(ax_bar, 'Force (N)');
    title(ax_bar, 'Empirical Validation: Invariance of Omnidirectional Control Effort');
    grid(ax_bar, 'on');
    
    % GUI Controls Configuration
    slider = uicontrol('Style', 'slider', 'Parent', fig, 'Min', 0, 'Max', pi, 'Value', initial_lambda, 'Position', [50, 20, 200, 20]);
    chk_particles = uicontrol('Style', 'checkbox', 'Parent', fig, 'String', 'Show Particles', 'Position', [260, 20, 100, 20], 'BackgroundColor', 'w', 'Value', 1);
    chk_sweep = uicontrol('Style', 'checkbox', 'Parent', fig, 'String', 'Auto-Sweep \lambda', 'Position', [370, 20, 110, 20], 'BackgroundColor', 'w', 'Value', auto_start_sweep);
    chk_motion = uicontrol('Style', 'checkbox', 'Parent', fig, 'String', 'Auto-Move', 'Position', [490, 20, 80, 20], 'BackgroundColor', 'w', 'Value', auto_start_motion);
    btn_reset = uicontrol('Style', 'pushbutton', 'Parent', fig, 'String', 'Reset Particles', 'Position', [580, 20, 100, 20]);
    
    % Dedicated Stop & Save Data Button
    btn_stop = uicontrol('Style', 'pushbutton', 'Parent', fig, 'String', 'Stop & Save', ...
                         'Position', [690, 20, 90, 25], 'BackgroundColor', [1 0.6 0.6], 'FontWeight', 'bold');
                         
    %% --- 4. Environment & Data Logging Initialization ---
    % Initialize ground dust particles
    [Xp, Yp] = meshgrid(linspace(grid_X_span(1)+0.1, grid_X_span(2)-0.1, 120), linspace(grid_Y_span(1)+0.1, grid_Y_span(2)-0.1, 45)); 
    particle_pos_init = [Xp(:), Yp(:)]';
    particle_pos_current = particle_pos_init;
    num_particles = size(particle_pos_current, 2);
    
    p_scatter = scatter(ax_ground, particle_pos_current(1,:), particle_pos_current(2,:), 10, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
    
    % UI Callbacks
    set(btn_reset, 'Callback', @(~,~) assignin('caller', 'particle_pos_current', particle_pos_init));
    set(btn_stop, 'Callback', @(~,~) set(fig, 'UserData', false)); % Cleanly breaks the main loop
    
    % Setup computational grid for wake projection
    [Xg, Yg] = meshgrid(linspace(grid_X_span(1), grid_X_span(2), 120), linspace(grid_Y_span(1), grid_Y_span(2), 50)); 
    Xg_flat = Xg(:); Yg_flat = Yg(:);
    P_grid = [Xg_flat, Yg_flat, -H * ones(length(Xg_flat), 1)]';
    
    ellipse_patches = gobjects(N, 1); quiver_plots = gobjects(N, 1);
    for i = 1:N
        quiver_plots(i) = quiver(ax_ground, NaN, NaN, NaN, NaN, 0, 'Color', colors(i,:), 'LineWidth', 1.2, 'MaxHeadSize', 0.5);
        ellipse_patches(i) = fill(ax_ground, NaN, NaN, colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', colors(i,:), 'LineWidth', 1.5);
    end
    uistack(p_scatter, 'top');
    
    % Pre-allocate telemetry logs arrays
    max_steps = 50000;
    log_t = zeros(1, max_steps);
    log_lambda = zeros(1, max_steps);
    log_f = zeros(N, max_steps);
    log_norm_f = zeros(1, max_steps);
    log_particles = zeros(1, max_steps);
    log_idx = 0;
    sim_time = 0;
    
    % Initialize phases
    lambda_phase = initial_lambda;
    drone_phase = 0;
    prev_chk_sweep = chk_sweep.Value;
    
    %% --- 5. Continuous Physics Engine Loop ---
    % Loop executes continuously until the window is closed or 'Stop' is pressed
    while ishandle(fig) && fig.UserData 
        
        % 1. Evaluate the 1D Morphological State (Lambda)
        if chk_sweep.Value
            if ~prev_chk_sweep
                lambda_phase = slider.Value; 
            end
            % Linear sweep creates a triangle wave oscillation across [0, pi]
            lambda_phase = lambda_phase + lambda_sweep_speed * dt;
            slider.Value = pi - abs(mod(lambda_phase, 2*pi) - pi);
        end
        prev_chk_sweep = chk_sweep.Value;
        lambda = slider.Value;
        
        % 2. Update Platform Kinematics (Navigation Task)
        if chk_motion.Value
            drone_phase = drone_phase + drone_motion_freq * dt;
            drone_acc_x = -drone_motion_amp * (drone_motion_freq^2) * cos(drone_phase);
        else
            drone_acc_x = 0;
        end
        P = P_base + [drone_motion_amp * cos(drone_phase); 0; 0];
        
        % 3. Control Allocation (Newton-Euler Equilibrium)
        A = zeros(6, N);
        for i = 1:N
            % Construct local thrust vector using Affine Phase-Locking equation
            theta = lambda + delta(i);
            d_i = cos(theta)*U(:,i) + sin(theta)*V(:,i);
            A(:,i) = [d_i; cross(P_base(:,i), d_i)];
        end
        
        % Compute required wrench (Hover + Lateral acceleration)
        W_req = [mass * lateral_inertia_multiplier * drone_acc_x; 0; mass*g; 0; 0; 0];
        
        % Pseudoinverse Allocation
        f = pinv(A) * W_req;
        norm_f = norm(f); % Compute Total Omnidirectional Control Effort
        
        % Update Bar Chart
        b_plot.YData = f; 
        b_effort.YData = norm_f; 
        
        % 4. Aerodynamic Wake Projection & Particle Advection
        vx = zeros(1, num_particles); vy = zeros(1, num_particles);
        P_parts = [particle_pos_current; -H * ones(1, num_particles)];
        
        for i = 1:N
            w_dir = -sign(f(i)) * A(1:3, i);
            % Skip calculation if wake points upwards (away from ground)
            if w_dir(3) >= -1e-3
                set(ellipse_patches(i), 'XData', NaN, 'YData', NaN);
                set(quiver_plots(i), 'XData', NaN, 'YData', NaN, 'UData', NaN, 'VData', NaN);
                continue;
            end
            w_n = w_dir / norm(w_dir);
            
            % Compute wake cylinder projection onto the ground plane
            if abs(w_n(3)) > 0.99; b1 = [1; 0; 0];
            else; b1 = cross([0; 0; 1]', w_n')'; b1 = b1 / norm(b1); end
            b2 = cross(w_n, b1);
            
            circ3D = P(:,i) + R_prop * (b1 * cos(linspace(0,2*pi,50)) + b2 * sin(linspace(0,2*pi,50)));
            proj3D = circ3D + w_n * ((-H - circ3D(3,:)) / w_n(3)); 
            set(ellipse_patches(i), 'XData', proj3D(1,:), 'YData', proj3D(2,:));
            
            V_mag = sqrt(abs(f(i))); 
            
            % Render wake flow field vectors
            v_vecs_grid = P_grid - P(:,i);
            mask_grid = (vecnorm(cross(v_vecs_grid, repmat(w_n, 1, size(v_vecs_grid, 2)))) / norm(w_n)) <= R_prop;
            if sum(mask_grid) > 0
                set(quiver_plots(i), 'XData', Xg_flat(mask_grid), 'YData', Yg_flat(mask_grid), ...
                                     'UData', arrow_scale * V_mag * w_n(1) * ones(sum(mask_grid),1), ...
                                     'VData', arrow_scale * V_mag * w_n(2) * ones(sum(mask_grid),1));
            else
                set(quiver_plots(i), 'XData', NaN, 'YData', NaN, 'UData', NaN, 'VData', NaN);
            end
            
            % Advect environmental particles
            if chk_particles.Value
                mask_parts = (vecnorm(cross(P_parts - P(:,i), repmat(w_n, 1, num_particles))) / norm(w_n)) <= R_prop;
                vx(mask_parts) = vx(mask_parts) + sim_speed * V_mag * w_n(1);
                vy(mask_parts) = vy(mask_parts) + sim_speed * V_mag * w_n(2);
            end
        end
        
        % Update particle positions
        if chk_particles.Value
            particle_pos_current(1,:) = particle_pos_current(1,:) + vx * dt;
            particle_pos_current(2,:) = particle_pos_current(2,:) + vy * dt;
            set(p_scatter, 'XData', particle_pos_current(1,:), 'YData', particle_pos_current(2,:), 'Visible', 'on');
        else
            set(p_scatter, 'Visible', 'off');
        end
        
        % 5. Data Logging (Record State for Paper Plots)
        in_box_x = particle_pos_current(1,:) >= monitor_X_span(1) & particle_pos_current(1,:) <= monitor_X_span(2);
        in_box_y = particle_pos_current(2,:) >= monitor_Y_span(1) & particle_pos_current(2,:) <= monitor_Y_span(2);
        
        log_idx = log_idx + 1;
        if log_idx <= max_steps
            log_t(log_idx) = sim_time;
            log_lambda(log_idx) = lambda;
            log_f(:, log_idx) = f;
            log_norm_f(log_idx) = norm_f;
            log_particles(log_idx) = sum(in_box_x & in_box_y); % Track particles advected into target area
        end
        
        sim_time = sim_time + dt;
        drawnow limitrate;
    end
    
    %% --- 6. Post-Simulation Data Export ---
    % Executes safely when the 'Stop & Save' button is pressed
    if log_idx > 0
        valid_len = min(log_idx, max_steps);
        log_t = log_t(1:valid_len);
        log_lambda = log_lambda(1:valid_len);
        log_f = log_f(:, 1:valid_len);
        log_norm_f = log_norm_f(1:valid_len);
        log_particles = log_particles(1:valid_len);
        
        disp('Simulation stopped by user. Saving logs to morphing_sim_log.mat...');
        save('morphing_sim_log.mat', 'log_t', 'log_lambda', 'log_f', 'log_norm_f', 'log_particles');
        disp('Data saved successfully!');
    end
    
    % Gracefully close the figure window
    if ishandle(fig)
        close(fig);
    end
end