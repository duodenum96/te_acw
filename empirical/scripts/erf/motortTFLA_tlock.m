%% Calculate motort erfs
% pid: 47004
clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

load Motort_subjs.mat
task = 'Motort';
loadext = 'TFLA';
taskext = 'TFLA_noprestim';
nsubj = size(subjs, 1);
scans = {'10', '11'};
nscan = length(scans);

badsubjs = [];
% Redefinition of trials
cfgx = [];
cfgx.toilim = [0 1050/1000];

for i = 1:nsubj
    isubj = subjs(i, :);
    try
        tic
        try
            cd(['/BICNAS2/ycatal/sourceparcellated/', isubj, '/', task, '/icamne/source'])
        catch
            cd(['/BICNAS2/ycatal/sourceparcellated/', isubj, '/', task, '/source'])
        end
        roidata = cell(1, nscan);
        cfg = [];
        for j = 1:nscan
            roidata{j} = load([isubj, '_MEG_', scans{j}, '-', task, '_', loadext, '_glasser.mat']);
            roidata{j} = ft_redefinetrial(cfgx, roidata{j}.roidata);
        end
        roidataboth = ft_appenddata(cfg, roidata{:});

        avgroidata = ft_timelockanalysis(cfg, roidataboth);
        save([isubj, '_MEG', '-', task, '_', taskext, '_erf.mat'], 'avgroidata')
        disp(['Subject ', num2str(i), ' finished!!!!'])
        toc
    catch
        disp(['Bad subj detected, id: ', subjs(i, :)])
        badsubjs = [badsubjs; subjs(i, :)];
    end
end
save(['/BICNAS2/ycatal/sourcereconstruction/data/', task, taskext, '_erf_badsubjs.mat'], 'badsubjs')
disp('Everything finished!!!!!!!!!! YAYYYYYYYYYY!!!!!!!!!!!!!')
