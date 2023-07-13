clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

load Restin_subjs.mat
nsubj = size(subjs, 1);
scans = {'3', '4', '5'};
nscan = length(scans);

task = 'Wrkmem';
taskext = 'TIM_noprestim';
loadext = 'TIM';

load(['/BICNAS2/ycatal/sourcereconstruction/data/ga/', task, taskext, '_selected_literature.mat'])
nroi = length(selected);

% Redefinition of trials
cfg = [];
cfg.length = 2;
cfg.overlap = 0;

badsubjs = [];
acw0 = zeros(nroi, nsubj, nscan);
acw50 = zeros(nroi, nsubj, nscan);

for i = 1:nsubj
    try
        tic
        try
            cd(['/BICNAS2/ycatal/sourceparcellated/', subjs(i, :), '/Restin/icamne/source'])
        catch
            cd(['/BICNAS2/ycatal/sourceparcellated/', subjs(i, :), '/Restin/source'])
        end
        for j = 1:nscan
            load([subjs(i, :), '_MEG_', scans{j}, '-Restin_glasser.mat'])
            cfgx = []; cfgx.continuous = 'yes';
            roidata = ft_redefinetrial(cfgx, roidata);
            roidata = ft_redefinetrial(cfg, roidata);

            data = cat(3, roidata.trial{:});
            data([1 182], :, :) = [];
            data = data(selected, :, :);
            
            [acw0(:, i, j), acw50(:, i, j)] = acw_3d(data, roidata.fsample);
        end
        disp(['Subject ', num2str(i), ' done.'])
        toc
    catch
        disp(['Bad subj detected, id: ', subjs(i, :)])
        badsubjs = [badsubjs; subjs(i, :)];
    end
end

disp('Everything finished!!!!!!!!!!!!!')
save(['/BICNAS2/ycatal/sourcereconstruction/data/acw_restmatched_', task, loadext, '_round2.mat'], 'acw50', 'acw0', 'badsubjs')