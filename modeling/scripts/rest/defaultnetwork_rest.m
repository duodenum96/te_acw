%% Do ACW-TE relationship with default parameters (VIS)
%% Start with resting state. 
% For task, I will do something like an ERP analysis using a boxcar input
cd /BICNAS2/ycatal/te_acw/modeling
init
irois = visual;

ntrials = 125; % Went forward in time and hardcoded to initialize qs. 
% This is true for p.dt = 1; p.tspan = 300;

acw0 = zeros(nvis, nsim);
acw50 = zeros(nvis, nsim);
acwdr = zeros(nvis, nsim);

fs = ((1/p.dt) * 1000) / 2; % *1000 for ms to second, / 2 for downsampling

for i = 1:nsim
    p.I_ext_E  = zeros([nroi ntime]);
    p.I_ext_E(1,:) = noises(i, :);
    p.I_ext_E(2:29, :) = squeeze(noisesrest(:, :, i));
    [v_E, time] = chaudhuri(p);
    v_E = downsample(v_E', 2)'; % Downsample to 500 Hz fs
    v_E = v_E(:, 1:(end-1)); % To make it divisible to 1000
    data = reshape(v_E, 29, 1000, []); % Divide into 2 second trials
    data = data(irois, :, :);
    assert(size(data, 3) == ntrials); % Double check if you change parameters
    [acw0(:, i), acw50(:, i)] = acw_3d(data, fs);
end
save('.gitignore/results/acws_defaultnetw.mat', 'acw0', 'acw50')
disp('Everything done!!!!!!!!!')
