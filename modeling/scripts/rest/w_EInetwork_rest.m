%% Do ACW-TE relationship with changing eta parameters
cd /BICNAS2/ycatal/te_acw/modeling
init
irois = visual;
nirois = nvis;
valuetoplay = 'w_EI';

ntrials = 125; % Went forward in time and initialized

c = [1 1.1 1.2];

acw0 = zeros(nirois, nsim, length(c));
acw50 = zeros(nirois, nsim, length(c));

fs = ((1/p.dt) * 1000) / 2; % *1000 for ms to second, / 2 for downsampling

defaultval = p.(valuetoplay);
savename = [valuetoplay, 'netw_rest'];

for ic = 1:length(c)
    p.(valuetoplay) = defaultval * c(ic);
    for i = 1:nsim
        tic
        p.I_ext_E  = zeros([nroi ntime]);
        p.I_ext_E(1,:) = noises(i, :);
        p.I_ext_E(2:29, :) = squeeze(noisesrest(:, :, i));
        [v_E, time] = chaudhuri(p);
    
        v_E = downsample((v_E'), 2)'; % Downsample to 500 Hz fs
        time = downsample(time, 2);
        v_E = v_E(:, 1:(end-1), :); % To make it divisible to 1000
        time = time(1:(end-1), :);

        data = reshape(v_E, 29, 1000, []);
        data = data(irois, :, :);
        
        assert(size(data, 3) == ntrials); % Double check if you change parameters

        [acw0(:, i, ic), acw50(:, i, ic)] = acw_3d(data, fs);
        disp(['Simulation ', num2str(i), ' done!!!!!!!'])
        toc
    end
    disp(['c ', num2str(ic), ' done!!!!!!!'])
end
save(['.gitignore/results/acws_', savename, '.mat'], 'acw0', 'acw50')
disp('Everything done!!!!!!!!!')
