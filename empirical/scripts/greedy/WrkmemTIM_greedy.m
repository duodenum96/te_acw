%% Do the greedy algorithm for each trials in each subject for each task in each trial definition

clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

task = 'Wrkmem';
taskext = 'TIM_noprestim';
loadext = 'TIM';
load([task, '_subjs.mat'])
load(['/BICNAS2/ycatal/sourcereconstruction/data/ga/', task, taskext, '_selected.mat'])
nsubj = size(subjs, 1);
scans = {'6', '7'};
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

            roidata = load([isubj, '_MEG-', task, '_', taskext, '_erf.mat']);
            roidata = roidata.avgroidata;
            ntime = length(roidata.time);
            datamatrix = roidata.avg;
            datamatrix([1 182], :) = [];

            datamatrix = datamatrix(selected, :);
            
            q_xx = zeros(nroi, nroi);
            q_xy = zeros(nroi, nroi);
            q_yx = zeros(nroi, nroi);
            q_yy = zeros(nroi, nroi);

	    disp(['Starting greedy for subject ', num2str(i), ...
                                ' scan ', scans{j}])
	    tic            
        for p = 1:nroi
            for q = p:nroi
                if p == q
                    q_xx(p, q) = 1; q_xy(p, q) = 1; 
                    q_yx(p, q) = 1; q_yy(p, q) = 1;
                else
                    x = datamatrix(p, :);
                    y = datamatrix(q, :);
                    [q_xx(p, q), q_xy(p, q), q_yx(p, q), q_yy(p, q)] = greedy(x, y, qmax);
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
