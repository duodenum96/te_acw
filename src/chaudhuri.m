function [v_E, time] = chaudhuri(p)
    %% Chaudhuri Model described in Chaudhuri et al. (2015)
    % Authored by Yasir Çatal a.k.a. Duodenum
    %% Initialize Free Parameters
    
    tau_E = p.tau_E;        % Intrinsic time constant(ms)
    tau_I = p.tau_I;
    beta_E = p.beta_E;      % Slope of firing rate (Hz / pA)
    beta_I = p.beta_I;
    w_EE = p.w_EE;          % Excitatory to excitatory coupling (pA / Hz)
    w_EI = p.w_EI;          % Inhibitory to excitatory
    w_IE = p.w_IE;
    w_II = p.w_II;
    mu_EE = p.mu_EE;        % Fixed parameter that control excitatory input (pA / Hz)
    mu_IE = p.mu_IE;
    eta = p.eta;            % Scaling parameter for hierarchy
    dt = p.dt;              % Time steps (ms)
    tspan = p.tspan;        % How many seconds
    h = p.h;                % Hierarchy thing
    J = p.J;                % SC matrix
    I_ext_E = p.I_ext_E;    % External current
    
    time = 0:dt:tspan/dt*1000;
    % Initialize Chaudhuri
    nroi = length(J);
    
    v_E = zeros([nroi length(time)]);                                              % Starting firing rate for excitatory, rows are regions, columns are time
    v_I = zeros([nroi length(time)]);
    zerovec = zeros([nroi 1]); % For the rectifier stuff
    
    %% Run the model
    for t = 1:length(time)-1
        I_lr_E = mu_EE .* sum(J' .* v_E(:,t))';
        I_lr_I = mu_IE .* sum(J' .* v_E(:, t))';
        I_E = (1 + eta.*h) .* (w_EE .* v_E(:, t) + I_lr_E) - w_EI * v_I(:, t) + I_ext_E(:, t);
        I_I = (1 + eta.*h) .* (w_IE * v_E(:, t) + I_lr_I) - w_II * v_I(:, t);
        v_E(:, t + 1) = v_E(:, t) + dt * (-v_E(:, t) + beta_E * max([I_E, zerovec], [], 2)) / tau_E;
        v_I(:, t + 1) = v_I(:, t) + dt * (-v_I(:, t) + beta_I * max([I_I, zerovec], [], 2)) / tau_I;
    end
    
    % Check for going to infinity
    if any(v_E(:) > 1)
        v_E = NaN([nroi length(time)]);
        time = NaN([nroi length(time)]);
        warning('Things are out of control, check the results carefully')
    end
    
    % Take the stable part - discard the first 50 seconds
    v_E = v_E(:, 50/dt*1000+1:end);
    time = time(50/dt*1000+1:end) / 1000; % Divide by 1000 to turn to seconds
    end