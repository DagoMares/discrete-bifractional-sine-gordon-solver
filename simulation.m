function simulation()
    % Parameters
    beta = 1.1; alpha = 1.1; gam = 0.5; lam = 1;%/(2*pi^2);
    a = -4; b = 4; % Grid limits
    tau = 0.01; h = 0.1; M = round(abs(b - a)/h); % Time-space steps
    nvar = (M + 1)^2; time=5;
    % Vid name
    %name = 'cos';
    name = 'soliton_t_0.005_h_0.05_a_1.5_b_1.5';
    % Frames 
    frames = 10;
    % Snapshots time
    sstimes = [0 5 10 15 20 25 30];
    % Snapshots and video z limits for numerical solution
    num_z_lim = [-3 3];
    % Initial conditions
    phi = zeros(nvar, 1);
    zk = zeros(nvar, 1);
    wk = zeros(nvar, 1);
    
    % Initialize wk and phi
    for j = M:-1:0
        for i = 0:M
            lin = (i + 1) + (M - j)*(M + 1);
            %wk(lin) = cos(pi*(a + i*h)) * cos(pi*(a + j*h));
            %wk(lin) = sin(2*pi*(a + i*h)) * sin(2*pi*(a + j*h));
            wk(lin) = 2*atan(exp(3-5*sqrt((a+i*h)^2+(a+j*h)^2)));
            phi(lin) = 1;
            % Threshold small values to zero
            if abs(wk(lin)) < 1e-10
                wk(lin) = 0;
            end
        end
    end
    
    % Define zero positions
    zeropos = [1:(M+1),(M+2):(M+1):(M+1)*(M+1), ...
               (2*(M+1)):(M+1):(M+1)*(M+1), ...
               (2 + (M+1)*M):(2 + (M+1)*M + M - 1)];
    wk(zeropos) = 0;
    zk(zeropos) = 0;
    
    % Precompute matrices
    LA = zeros(nvar); LB = zeros(nvar);
    Lmatrix = zeros(nvar); Rmatrix = zeros(nvar);
    I = eye(M + 1);
    r = (2 + gam*tau)/tau^2; R = r * I;
    
    % Compute H matrices
    Halpha = H_m(M, alpha, h);
    Hbeta_row = H_m(M, beta, h); Hbeta = Hbeta_row(1, :);
    Ha2 = H_m(M, alpha/2, h); 
    Hb2_row = H_m(M, beta/2, h); Hb2 = Hb2_row(1, :);
    
    % Construct D matrices
    D1 = Halpha + Hbeta(1)*I - (2/lam)*R;
    D2 = Halpha + Hbeta(1)*I + (2/lam)*R;
    DB = Hb2(1)*I;
    
    % Build L1I, L1D, L1B blocks
    L1I = D1/2; L1D = D2/2; L1B = DB/2;
    for i = 2:(M+1)
        L1I = [L1I, Hbeta(i)*I]; %#ok<AGROW>
        L1D = [L1D, Hbeta(i)*I]; %#ok<AGROW>
        L1B = [L1B, Hb2(i)*I]; %#ok<AGROW>
    end
    
    % Assemble global matrices
    for i = 1:(M+1)
        group = 1 + (i-1)*(M+1);
        cols = group:((M+1)^2);
        slice_cols = 1:((M+2 - i)*(M+1));
        Lmatrix(group:group+M, cols) = L1I(:, slice_cols);
        Rmatrix(group:group+M, cols) = L1D(:, slice_cols);
        LA(group:group+M, group:group+M) = Ha2;
        LB(group:group+M, cols) = L1B(:, slice_cols);
    end
    
    % Finalize matrices
    Lmatrix = (-lam/2) * (Lmatrix + Lmatrix');
    Rmatrix = (lam/2) * (Rmatrix + Rmatrix');
    LB = LB + LB';
    Linvmatrix = inv(Lmatrix);
    
    % Time stepping setup
    %time = 1;
    timesteps = ceil(time/tau);
    E = zeros(timesteps, 1);
    evolution = cell(timesteps, 1);
    
    % Main simulation loop
    for l = 1:timesteps
        % Energy computation
        da = sum((LA * wk).^2);
        db = sum((LB * wk).^2);
        vk2 = sum(zk.^2);
        E(l) = 0.5 * h^2 * (vk2 + lam*(da + db)) + h^2 * sum(phi .* (1 - cos(wk)));
        
        % Reshape wk into spatial matrix
        U = reshape(wk, M+1, M+1)';
        evolution{l} = U;
        
        % Compute forcing term
        %Fk = Fkfunction(l-1, nvar, h, tau, M, a);
        Fk = 0;
        % Nonlinear iteration
        wkm = 2 * ones(size(wk));
        err = 1; tol = 1e-8; max_iter = 150; counter = 1;
        while err > tol && counter <= max_iter
            delta = wkm - wk;
            delta(delta == 0) = eps; % Avoid division by zero
            wkm1 = Linvmatrix * (Rmatrix * wk + (2/tau) * zk + ...
                   phi .* (cos(wkm) - cos(wk)) ./ delta + Fk); %#ok<MINV>
            err = sum(abs(wkm1.^2 - wkm.^2));
            wkm = wkm1;
            counter = counter + 1;
        end
        wkm1(zeropos) = 0; 
        zk = 2*(wkm1 - wk)/tau - zk;
        zk(zeropos) = 0;
        wk = wkm1;
    end
    
    % Plot results
    x = a:h:b;
    y = b:-h:a;
    [X, Y] = meshgrid(x, y);
    
    % Create animation with comparison
    %create_comparison_animation(X, Y, evolution, tau, a, b, h,name, frames, num_z_lim);
    
    % Calculate and plot norm differences
    %plot_norm_differences(X, Y, evolution, tau, timesteps,h);

    % Plot energy evolution
    plot_energy_solution(E, tau, timesteps);

    % Solution snapshots
    %save_solution_snapshots(evolution, X, Y, tau, a, b, h, sstimes, alpha, beta, num_z_lim);

    % Unique animation
    %animate_pde_solution(evolution, x, y, tau, name, frames, num_z_lim);

end

% New function for energy visualization
function plot_energy_solution(E, tau, timesteps)
    % Create time vector
    time_vector = (0:tau:(timesteps-1)*tau);
    % Mean energy
    mE = mean(E);
    % Maximum energy - minimum energy
    maxE = max(E);     minE = min(E);
    % Create energy figure
    figure('Position', [100 100 800 400]);
    
    % ===== MODIFY THIS SECTION =====
    % Original plot command:
    plot(time_vector, E, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980],...
        'DisplayName','Instantaneous energy');
    
    % Example 1: Change to blue dashed line
    % plot(time_vector, E, '--', 'LineWidth', 2, 'Color', 'blue');
    
    % Example 2: Add markers
    % plot(time_vector, E, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', '#0072BD');
    
    % Example 3: Add moving average (uncomment to use)
    % window_size = 20;
    % mov_avg = movmean(E, window_size);
    % hold on;
    % plot(time_vector, mov_avg, 'k-', 'LineWidth', 2);
    % hold off;
    % ==============================
    hold on;
    plot(xlim, [mE mE], 'r--', 'LineWidth', 1.5,...
        'DisplayName', sprintf(' Mean energy = %.4f', mE)); % Mean line
    plot(xlim, [max(E) max(E)], 'b--', 'LineWidth', 1.5,...
        'DisplayName', sprintf(' Minimal energy = %.4f', minE)); % Max line
    plot(xlim, [min(E) min(E)], 'b--', 'LineWidth', 1.5,...
        'DisplayName', sprintf(' Maximal energy = %.4f', maxE)); % Min line
    legend('Location', 'best', 'FontSize', 10);
    hold off;

    xlabel('Time (t)', 'FontSize', 12);
    ylabel('Energy', 'FontSize', 12);
    title('Energy evolution of the solution', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12);
    
    
end

% Function for norm difference visualization
function plot_norm_differences(X, Y, evolution, tau, timesteps,h)
    % Calculate time vector
    time_vector = (0:tau:(timesteps-1)*tau);
    
    % Preallocate norm array
    norms = zeros(timesteps, 1);
    
    % Calculate norms at each timestep
    for l = 1:timesteps
        t = (l-1)*tau;
        U_pde = evolution{l};
        U_analytical = cos(pi*X) .* cos(pi*Y) .* cos(t);
        su = sum(sum((U_pde - U_analytical).^2));
        norms(l) = sqrt(h^2*su);%norm(U_pde - U_analytical, 'fro'); % Frobenius norm
    end
    
    % Create norm difference figure
    figure('Position', [100 100 800 400]);
    plot(time_vector, norms, 'LineWidth', 2);
    xlabel('Time (t)', 'FontSize', 12);
    ylabel('Discrete L-2 error', 'FontSize', 12);
    title('Numerical vs Analytical solution', 'FontSize', 14);
    grid on;
    set(gca, 'FontSize', 12);
    
    % Add reference line for scale
    hold on;
    plot(xlim, [0 0], 'k--');
    hold off;
end

% Animation function with analytical comparison
function create_comparison_animation(X, Y, evolution, tau, a, b, ~,name,frames, num_z_lim)
    %Set up video writer
    video_name = name;
    writer = VideoWriter(video_name, 'MPEG-4');
    writer.FrameRate = frames;
    open(writer);
    
    % Create figure with 2 subplots
    fig = figure('Color', 'white', 'Position', [100 100 1200 600]);
    
    % Fixed axis settings ️
    z_limits = num_z_lim;  % Fixed z-axis range for numerical solution
    ana_z_limits = [-1 1]; % Fixed z-axis range for analytical solution
    xy_limits = [a b];  % Fixed x/y-axis range
    
    % Custom color palette
    colormap(fig, "cool"); % Change to preferred palette
    
    % Calculate analytical solution bounds
    %t_max = (length(evolution)-1)*tau;
    [~, ~] = meshgrid(linspace(a, b, size(X,2)), linspace(a, b, size(Y,1))');
    
    % Set consistent color limits
    all_pde = cellfun(@(x) x(:), evolution, 'UniformOutput', false);
    all_pde = vertcat(all_pde{:});
    c_limits = [min(all_pde), max(all_pde)];

    for k = 1:length(evolution)
        % Current time
        t = (k-1)*tau;
        
        % PDE Solution plot
        subplot(1,2,1);
        surf(X, Y, evolution{k});
        title(sprintf('PDE Solution\ntime = %.2f', t));
        axis([xy_limits xy_limits z_limits]); % Fixed 3D axes ️
        clim(c_limits);
        xlabel('x'); ylabel('y'); zlabel('U');
        shading interp;
        view(-30, 5);
        
        % Analytical solution plot
        subplot(1,2,2);
        U_analytical = cos(pi*X) .* cos(pi*Y) .* cos(t);
        %U_analytical = 4*atan(sinh(X/2+Y/2).*sin(t/sqrt(2))) ;
        surf(X, Y, U_analytical);
        title(sprintf('Analytical solution\nu(x,y,t) = cos(πx)cos(πy)cos(t)'));
        %title(sprintf('Analytical Solution'));
        axis([xy_limits xy_limits ana_z_limits]); % Fixed 3D axes ️
        clim(c_limits);
        xlabel('x'); ylabel('y'); zlabel('U');
        shading interp;
        view(-30, 5);
        
        % Capture frame
        drawnow;
        frame = getframe(fig);
        writeVideo(writer, frame);
    end
    
    close(writer);
    fprintf('Comparison animation saved as %s\n', video_name);
end

% Helper function for Fk
function Ff = Fkfunction(index, nvar, h, tau, M, a)
    Ff = zeros(nvar, 1);
    for j = M:-1:0
        for i = 0:M
            lin = (i + 1) + (M - j)*(M + 1);
            term = cos(pi*(a + i*h)) * cos(pi*(a + j*h)) * cos((index+0.5) * tau);
            %term = sin(2*pi*(a + i*h)) * sin(2*pi*(a + j*h)) * cos(sqrt(8)*pi*(index+0.5) * tau);
            Fk = sin(term);
            %term_next = cos(pi*(a + i*h)) * cos(pi*(a + j*h)) * cos((index + 1) * tau);
            %term_next = sin(2*pi*(a + i*h)) * sin(2*pi*(a + j*h)) * cos(sqrt(8)*pi*(index + 1) * tau);
            %Fk1 = sin(term_next);
            Ff(lin) = (Fk); %+ Fk1)/2;
            if abs(Ff(lin)) < 1e-10
                Ff(lin) = 0;
            end
        end
    end
end

% Gamma coefficient function
function gal = g_l_a(alpha, l)
    gal = zeros(1, l+1);
    gal(1) = gamma(alpha + 1) / (gamma(alpha/2 + 1)^2);
    for k = 2:(l+1)
        gal(k) = (1 - (alpha + 1)/(alpha/2 + (k-2) + 1)) * gal(k-1);
    end
end

% H matrix generator
function H = H_m(M, alpha, h)
    H = zeros(M+1);
    v = g_l_a(alpha, M);
    v(1) = v(1)/2;
    for i = 1:(M+1)
        H(i, i:end) = v(1:end - i + 1);
    end
    H = (H + H') / (-h^alpha);
end

% Snapshots
function save_solution_snapshots(evolution, X, Y, tau, a, b, ~, sstimes, alpha, beta, num_z_lim)
    % Times to capture (must be within simulation range)
    snapshot_times = sstimes;%[0.8, 1.5, 3.0];
    
    % Create analytical grid
    [~, ~] = meshgrid(linspace(a, b, size(X,2)), linspace(a, b, size(Y,1))');
    
    % Set consistent color limits
    all_pde = cellfun(@(x) x(:), evolution, 'UniformOutput', false);
    all_pde = vertcat(all_pde{:});
    c_limits = [min(all_pde), max(all_pde)];
    
    for t_target = snapshot_times
        % Calculate time index
        time_index = round(t_target/tau) + 1;
        
        % Skip if out of range
        if time_index > length(evolution)
            warning('Time %.1f is beyond simulation range. Skipping.', t_target);
            continue;
        end
        
        % Get numerical solution
        U_pde = evolution{time_index};
        
        % Compute analytical solution
        U_analytical = cos(pi*X) .* cos(pi*Y) .* cos(t_target);
        
        %% Create and save numerical solution figure
        fig_num = figure('Visible', 'off', 'Color', 'white', 'Position', [100 100 800 600]);
        surf(X, Y, U_pde);
        title(sprintf('Numerical solution')); %at t = %.1f', t_target));
        ax = gca;
        ax.TitleFontSizeMultiplier = 2;
        ax.FontSize = 8;
        zlim(num_z_lim);
        clim(c_limits);
        shading interp;
        view(-30, 15);
        xlabel('x'); ylabel('y'); zlabel('U');
        colormap('jet');
        colorbar;
        set(gca, 'FontSize', 12);
        
        % Save as EPS
        saveas(fig_num, sprintf('n_sol_t_%.1f_a_%.1f_b_%.1f.eps', t_target, alpha, beta), 'epsc');
        close(fig_num);
        
        %% Create and save analytical solution figure
        %{
        fig_ana = figure('Visible', 'off', 'Color', 'white', 'Position', [100 100 800 600]);
        surf(X, Y, U_analytical);
        title(sprintf('Analytical solution'));% at t = %.1f', t_target));
        ax = gca;
        ax.TitleFontSizeMultiplier = 2;
        ax.FontSize = 8;
        zlim([-1 1]);
        clim(c_limits);
        shading interp;
        view(-30, 15);
        xlabel('x'); ylabel('y'); zlabel('U');
        colormap('jet');
        colorbar;
        set(gca, 'FontSize', 12);
        
        % Save as EPS
        saveas(fig_ana, sprintf('a_sol_t_%.1f.eps', t_target), 'epsc');
        close(fig_ana);
        %}
    end
end

function animate_pde_solution(evolution, x, y, tau, name, frames, num_z_lim)
    % Create figure
    fig = figure('Color', 'white', 'Position', [100 100 800 600]);
    
    % Set fixed z-axis limits
    z_limits = num_z_lim;
    
    % Set custom color palette
    colormap(fig, 'jet');
    
    % Create initial surface plot
    [X, Y] = meshgrid(x, y);
    s = surf(X, Y, evolution{1});
    title(sprintf('Time = %.2f', 0));
    zlim(z_limits);
    shading interp;
    colorbar;
    xlabel('x'); ylabel('y'); zlabel('U');
    view(-30, 10);
    
    % Set up video writer
    video_name = name;%'pde_solution_animation.mp4';
    writer = VideoWriter(video_name, 'MPEG-4');
    writer.FrameRate = frames;
    open(writer);
    
    % Capture first frame
    drawnow;
    frame = getframe(fig);
    writeVideo(writer, frame);
    
    % Animate through time steps
    for k = 1:length(evolution)
        % Update surface data
        s.ZData = evolution{k};
        
        % Update title
        title(sprintf('Time = %.2f', (k-1)*tau));
        
        % Update plot
        drawnow;
        
        % Capture frame
        frame = getframe(fig);
        writeVideo(writer, frame);
    end
    
    close(writer);
    fprintf('Solution animation saved as %s\n', video_name);
end