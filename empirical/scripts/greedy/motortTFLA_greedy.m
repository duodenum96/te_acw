%% Do the greedy algorithm for each trials in each subject for each task in each trial definition

clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

task = 'Motort';
taskext = 'TFLA_noprestim';
loadext = 'TFLA';
load([task, '_subjs.mat'])
load(['/BICNAS2/ycatal/sourcereconstruction/data/ga/', task, taskext, '_selected.mat'])
nsubj = size(subjs, 1);
scans = {'10', '11'};
nscan = length(scans);

badsubjs = [];
cfg = [];
nroi = length(selected); % 10

% Redefinition of trials
cfgx = [];
cfgx.toilim = [0 2];

qmax = 20;

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

        for j = 1:nscan
            roidata = load([isubj, '_MEG_', scans{j}, '-', task, '_', loadext, '_glasser.mat']);
            roidata = ft_redefinetrial(cfgx, roidata.roidata);
            ntrials = length(roidata.trial);
            ntime = size(roidata.trial{1}, 2);
            datamatrix = zeros(nroi, ntime, ntrials);

            for k = 1:ntrials
                tmp = roidata.trial{k};
                tmp([1 182], :) = [];
                datamatrix(:, :, k) = tmp(selected, :);
            end
            
            q_xx = zeros(nroi, nroi, ntrials);
            q_xy = zeros(nroi, nroi, ntrials);
            q_yx = zeros(nroi, nroi, ntrials);
            q_yy = zeros(nroi, nroi, ntrials);

	    disp(['Starting greedy for subject ', num2str(i), ...
                  ' scan ', scans{j}])
            tic
            for k = 1:ntrials
                for p = 1:nroi
                    for q = p:nroi
                        if p == q
                            q_xx(p, q, k) = 1; q_xy(p, q, k) = 1; 
                            q_yx(p, q, k) = 1; q_yy(p, q, k) = 1;
                        else
                            x = datamatrix(p, :, k);
                            y = datamatrix(q, :, k);
                            [q_xx(p, q, k), q_xy(p, q, k), q_yx(p, q, k), q_yy(p, q, k)] = greedy(x, y, qmax);
                        end
                    end
                end
            end
	    toc
            savefile = [isubj, '_MEG', '-', task, '_', scans{j}, '_', taskext, '_qs.mat'];
            save(savefile, 'datamatrix', 'q_xx', 'q_xy', 'q_yx', 'q_yy')
            disp(['Subject ', num2str(i), ' done!!!!!!!!!!!!!!!!!'])
        end
    catch
        disp(['Bad subj detected, id: ', subjs(i, :)])
        badsubjs = [badsubjs; subjs(i, :)];
    end
end
