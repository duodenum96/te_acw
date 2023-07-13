%% Do source reconstruction for rest
clear, clc, close all
cd /BICNAS2/ycatal/sourcereconstruction
set_paths
%% Loop over subjects and do source reconstruction
load Restin_subjs.mat
load([projectfolder, '/data/badsubjs_Restin.mat'])


for subject = 1:size(subjs, 1)
    currentsubject = subject; % Subject

    subjectid = strtrim(subjs(currentsubject, :));
    scans = ['3', '4', '5'];
    task = 'Restin';
    experimentid = 'MEG';
    
    try
        for ii = 1:length(scans)
            % Clear stuff
            clear atlas c cc cfg channels comp_class datacat grad gradBalanced ...
                gridLF headmodel i inside_indices k mixing noise_level options roidata ...
                source sourcemodel2d sourcemodelsubj voxeldata w u s v
        
            % Initialize
            scanid = scans(ii);
            disp(['!!!!!!!!!!!!!!!!!!!!!!!!!Starting subject ', subjectid, ' Scan ', scanid])
            pipelinedatadir = ['/HCP/MEG/megdata/', subjectid, '/', experimentid];
            filename = ['/HCP/MEG/megdata/', subjectid, '/unprocessed/MEG/', scanid, '-', task, '/4D/c,rfDC'];
            outfolder = ['/BICNAS2/ycatal/sourceparcellated/', subjectid, '/', task];
            if ~exist(outfolder, 'dir')
                mkdir(outfolder)
            end
            % Run scripts
            hcp_rmegpreproc
            hcp_icamne
            disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!icamne done, scanid: ', num2str(scanid)])
            sourceparcellatev2
            disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!sourceparcellate done, scanid: ', num2str(scanid)])
        end
    catch
        disp(['!!!!!!!!!!!!!!!!!!bad subject detected, id: ', subjectid])
        badsubjs_Restin = [badsubjs_Restin; string(subjectid)];
        save([projectfolder, '/data/badsubjs_Restin.mat'], 'badsubjs_Restin')
    end
    disp(['Subject ', num2str(currentsubject), ' done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
end

disp('EVERYTHING FINISHED!!!')

% Subject 89 done
%% Go over

% diary badsubjs8_32.log
% % find(strtrim(string(subjs)) == "156334")
% badsubjs2 = [];
% for subject = 8:32
%     currentsubject = subject; % Subject
    
% %     subjectid = strtrim(subjs(currentsubject, :));
%     subjectid = char(badsubjs(subject));
%     scans = ['3', '4', '5'];
%     task = 'Restin';
%     experimentid = 'MEG';
%     try
%         for ii = 1:length(scans)
%             % Clear stuff
%             clear atlas c cc cfg channels comp_class datacat grad gradBalanced ...
%                 gridLF headmodel i inside_indices k mixing noise_level options roidata ...
%                 source sourcemodel2d sourcemodelsubj voxeldata w u s v
    
%             % Initialize
%             scanid = scans(ii);
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!Starting subject ', subjectid, ' Scan ', scanid])
%             pipelinedatadir = ['D:\HCPMEGData\', subjectid, '\', experimentid];
%             filename = ['D:\HCPMEGData\', subjectid, '\unprocessed\MEG\', scanid, '-', task, '\4D\c,rfDC'];
%             % Run scripts
%             hcp_rmegpreproc
%             hcp_icamne
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!icamne done, scanid: ', num2str(scanid)])
%             sourceparcellate
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!sourceparcellate done, scanid: ', num2str(scanid)])
%         end
%         disp(['Subject ', subjectid, ' DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%     catch
%         badsubjs2 = [badsubjs2, string(subjectid)];
%     end
% end
% diary off
% %% Check bad subjects bc I was stupid

% badsubjs2 = [];
% for subject = 1:length(badsubjs)
%     for scan = 5
%         if ~isfile("D:\HCPMEGData\" + badsubjs(subject) + "\MEG\Restin\icamne\" + badsubjs(subject) + "_MEG_" + scan + "-Restin_glasser.mat")
%             badsubjs2 = [badsubjs2, badsubjs(subject)];
%         end
%     end
% end
% %% Run badsubjs2 after redownloading raw data

% diary badsubjslast7.log
% % find(strtrim(string(subjs)) == "156334")
% badsubjs3 = [];
% for subject = 1:length(badsubjs2) % Start from 5 tomorrow
% %     currentsubject = subject; % Subject
    
% %     subjectid = strtrim(subjs(currentsubject, :));
%     subjectid = char(badsubjs2(subject));
%     scans = ['3', '4', '5'];
%     task = 'Restin';
%     experimentid = 'MEG';
%     try
%         for ii = 1:length(scans)
%             % Clear stuff
%             clear atlas c cc cfg channels comp_class datacat grad gradBalanced ...
%                 gridLF headmodel i inside_indices k mixing noise_level options roidata ...
%                 source sourcemodel2d sourcemodelsubj voxeldata w u s v
    
%             % Initialize
%             scanid = scans(ii);
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!Starting subject ', subjectid, ' Scan ', scanid])
%             pipelinedatadir = ['D:\HCPMEGData\', subjectid, '\', experimentid];
%             filename = ['D:\HCPMEGData\', subjectid, '\unprocessed\MEG\', scanid, '-', task, '\4D\c,rfDC'];
%             % Run scripts
%             hcp_rmegpreproc
%             hcp_icamne
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!icamne done, scanid: ', num2str(scanid)])
%             sourceparcellate
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!sourceparcellate done, scanid: ', num2str(scanid)])
%         end
%         disp(['Subject ', subjectid, ' DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
%     catch
%         disp(['Bad subject: ', subjectid])
%         badsubjs3 = [badsubjs3, string(subjectid)];
%     end
% end
% diary off
% %% Final loop
% subjectid = '181232';
%     task = 'Restin';
%     experimentid = 'MEG';
%             clear atlas c cc cfg channels comp_class datacat grad gradBalanced ...
%                 gridLF headmodel i inside_indices k mixing noise_level options roidata ...
%                 source sourcemodel2d sourcemodelsubj voxeldata w u s v
    
%             % Initialize
%             scanid = '5';
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!Starting subject ', subjectid, ' Scan ', scanid])
%             pipelinedatadir = ['D:\HCPMEGData\', subjectid, '\', experimentid];
%             filename = ['D:\HCPMEGData\', subjectid, '\unprocessed\MEG\', scanid, '-', task, '\4D\c,rfDC'];
%             % Run scripts
%             hcp_rmegpreproc
%             hcp_icamne
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!icamne done, scanid: ', num2str(scanid)])
%             sourceparcellate
%             disp(['!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!sourceparcellate done, scanid: ', num2str(scanid)])