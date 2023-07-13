%% Do source reconstruction for working memory task
clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths
%% Loop over subjects and do source reconstruction
load Motort_subjs.mat
% load([projectfolder, '/data/badsubjs_Wrkmem.mat'])
badsubjs_Motort = [];
scans = {'10', '11'};
task = 'Motort';
experimentid = 'MEG';

for subject = 1:size(subjs, 1)
    currentsubject = subject; % Subject

    subjectid = strtrim(subjs(currentsubject, :));
    
    
    try
        for ii = 1:length(scans)
            % Clear stuff
            clear atlas c cc cfg channels comp_class datacat grad gradBalanced ...
                gridLF headmodel i inside_indices k mixing noise_level options roidata ...
                source sourcemodel2d sourcemodelsubj voxeldata w u s v
        
            % Initialize
            scanid = scans{ii};
            disp(['!!!!!!!!!!!!!!!!!!!!!!!!!Starting subject ', subjectid, ' Scan ', scanid])
            pipelinedatadir = ['/HCP/MEG/megdata/', subjectid, '/', experimentid];
            filename = ['/HCP/MEG/megdata/', subjectid, '/unprocessed/MEG/', scanid, '-', task, '/4D/c,rfDC'];
            outfolder = ['/BICNAS2/ycatal/sourceparcellated/', subjectid, '/', task];

            if ~exist(outfolder, 'dir')
                mkdir(outfolder)
            end
            % Run scripts
            hcp_tmegpreproc
            hcp_icamne
            disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!icamne done, scanid: ', num2str(scanid)])
            for jj = 1:length(allTrlCfgs)
                taskextension = allTrlCfgs{jj}.trialdef.lockMnem;
                sourceparcellatev2
                disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!sourceparcellate done, scanid: ', num2str(scanid)])
            end
        end
    catch
        disp(['!!!!!!!!!!!!!!!!!!bad subject detected, id: ', subjectid])
        badsubjs_Motort = [badsubjs_Motort; string(subjectid)];
        save([projectfolder, '/data/badsubjs_Motort.mat'], 'badsubjs_Motort')
    end
    disp(['Subject ', num2str(currentsubject), ' done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
end

