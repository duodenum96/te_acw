%% Grand average of task - task extension
clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

load Wrkmem_subjs.mat
task = 'Wrkmem';
taskext = 'TIM_noprestim';
nsubj = size(subjs, 1);
scans = {'6', '7'};
nscan = length(scans);

badsubjs = [];
cfg = [];

alldata = cell(1, nsubj);
for i = 1:nsubj
    isubj = subjs(i, :);
    try
        tic
        try
            cd(['/BICNAS2/ycatal/sourceparcellated/', isubj, '/', task, '/icamne/source'])
        catch
            cd(['/BICNAS2/ycatal/sourceparcellated/', isubj, '/', task, '/source'])
        end
        
        data = load([isubj, '_MEG', '-', task, '_', taskext, '_erf.mat']);
        data.cfg = [];
        alldata{i} = data.avgroidata;
        disp(['Subject ', num2str(i), ' finished!!!!'])
        toc
    catch
        disp(['Bad subj detected, id: ', subjs(i, :)])
        badsubjs = [badsubjs; subjs(i, :)];
    end
end
granderf = ft_timelockgrandaverage(cfg, alldata{:});
granderf.cfg = [];
clear alldata

save(['/BICNAS2/ycatal/sourcereconstruction/data/ga/', task, taskext, '_ga.mat'], 'badsubjs', 'granderf', '-v7.3')
disp('Everything finished!!!!!!!!!! YAYYYYYYYYYY!!!!!!!!!!!!!')