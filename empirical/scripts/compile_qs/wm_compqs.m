%% Compile qs
clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

task = 'Wrkmem';
taskext = 'TIM_noprestim';
loadext = 'TIM';
load([task, '_subjs.mat'])
nsubj = size(subjs, 1);

load(['/BICNAS2/ycatal/sourcereconstruction/data/ga/', task, taskext, '_selected_literature.mat'])
nroi = length(selected); % 10

q_xx = zeros(nroi, nroi, nsubj);
q_xy = zeros(nroi, nroi, nsubj);
q_yx = zeros(nroi, nroi, nsubj);
q_yy = zeros(nroi, nroi, nsubj);

badsubjs = [];

for i = 1:nsubj
    isubj = subjs(i, :);
    try
        tic
        try
            savefolder = ['/BICNAS2/ycatal/sourceparcellated/', isubj, '/', task, '/icamne/source/'];
            cd(savefolder)
        catch
            savefolder = ['/BICNAS2/ycatal/sourceparcellated/', isubj, '/', task, '/source/'];
            cd(savefolder)
        end
            tedata = load([isubj, '_MEG-', task, '_', taskext, '_qs_erfall_lit.mat']);
            q_xx(:, :, i) = tedata.q_xx;
            q_xy(:, :, i) = tedata.q_xy;
            q_yx(:, :, i) = tedata.q_yx;
            q_yy(:, :, i) = tedata.q_yy;
    catch
        disp(['Bad subj detected, id: ', subjs(i, :)])
        badsubjs = [badsubjs; subjs(i, :)];
    end
end
savefile = ['/BICNAS2/ycatal/sourcereconstruction/data/qs/qs_', task, '.mat'];
save(savefile, 'q_xx', 'q_xy', 'q_yx', 'q_yy')
disp('Everything done ;)')