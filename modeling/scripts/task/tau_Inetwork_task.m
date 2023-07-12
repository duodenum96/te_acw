%% Do ACW-TE relationship with changing eta parameters
cd /BICNAS2/ycatal/te_acw/modeling
init
irois = visual;
nirois = nvis;
valuetoplay = 'tau_I';

starttime = 10;
dur = 2;
iti = dur + 1;
amplitude = 0.5;
dt = p.dt;
endtime = p.tspan;
% ntime
[box, offsets] = get_box(starttime, iti, dur, amplitude, endtime, dt, ntime);
ntrials = 83; % Hardcoded and checked in line 33

c = [0.6 0.8 1 1.2 1.4 1.6 1.8 2];

acw0 = zeros(nirois, nsim);
acw50 = zeros(nirois, nsim);
acwdr = zeros(nirois, nsim);
q_xx = zeros(nirois, nirois, nsim);
q_xy = zeros(nirois, nirois, nsim);
q_yx = zeros(nirois, nirois, nsim);
q_yy = zeros(nirois, nirois, nsim);
fs = ((1/p.dt) * 1000) / 2; % *1000 for ms to second, / 2 for downsampling

defaultval = p.(valuetoplay);
savename = [valuetoplay, 'netw_task'];

for ic = 1:length(c)
    p.(valuetoplay) = defaultval * c(ic);
    for i = 1:nsim
        tic
        p.I_ext_E  = zeros([nroi ntime]);
        p.I_ext_E(1,:) = noises(i, :) + box;
        p.I_ext_E(2:29, :) = squeeze(noisesrest(:, :, i));
        [v_E, time] = chaudhuri(p);
        [v_E, time] = trialize(v_E, time, offsets, dur, dt);
    
        v_E = pagetranspose(downsample(pagetranspose(v_E), 2)); % Downsample to 500 Hz fs
        time = downsample(time, 2);
        v_E = v_E(:, 1:(end-1), :); % To make it divisible to 1000
        time = time(1:(end-1), :);
        data = v_E(irois, :, :);
        assert(size(data, 3) == ntrials); % Double check if you change parameters
        [acw0(:, i, ic), acw50(:, i, ic)] = acw_3d(data, fs);
        % acwdr(:, i) = acw_dr_3d(data, fs, true, 1);
        
        datamatrix(:, :, i, ic) = mean(data, 3); % ERPs
    
        
        for j = 1:nirois
            for q = 1:nirois
                if j == q
                    q_xx(j, q, i, ic) = 1; q_xy(j, q, i, ic) = 1;
                    q_yx(j, q, i, ic) = 1; q_yy(j, q, i, ic) = 1;
                else
                    x = datamatrix(j, :, i, ic);
                    y = datamatrix(q, :, i, ic);
                    [q_xx(j, q, i, ic), q_xy(j, q, i, ic), ...
                        q_yx(j, q, i, ic), q_yy(j, q, i, ic)] = greedy(x, y, qmax);
                end
            end
        end
        disp(['Simulation ', num2str(i), ' done!!!!!!!'])
        toc
    end
    disp(['c ', num2str(ic), ' done!!!!!!!'])
end
save(['.gitignore/results/greedy_', savename, '.mat'], 'datamatrix', 'q_xx', 'q_xy', 'q_yx', 'q_yy')
save(['.gitignore/results/acws_', savename, '.mat'], 'acw0', 'acw50')
disp('Everything done!!!!!!!!!')
