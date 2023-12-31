function [cleanTrl] = hcp_extract_allfromrun2(inCfg)
%% This function performs the core processing of the tmegpreproc and rmegpreproc pipelines.
% It extracts  trial data for a give data group and fuses them with the
% results from hcp_baddata.m pipeline so that bad channels and trials that
% coincide with noisy periods are removed. In the case that a trial spans
% a long block of data, in order to avoid removing it entirely, the bad
% segments in the trials are replaced with nan. (This is the case for the Story/Math data groups BSENT and BUN).
% Then the IC components identified as related to heart or eye activity by hcp_icaclass.m pipeline are removed.
% Then the data is resampled to one fourth of the original sampling frequency in order to reduce the size of the
% dataset.
%
%
% INPUT:
%-------------------------------------------------------------------
% inCfg : This is a structure containing required parameters for the
%          analysis
%          Fields:
%                .datafile:  This is the raw data filename for a given scan.
%                .trl:       This is the trials definition matrix. It has 3 columns and number of rows equal to number of trials. The
%                              first column is the start sample, the second is the end sample and time offset of the start of the trial
%                              relative to the 0 reference point. This information is created by the trial definition functions.
%                .badchanfile: This is the file containing information about bad channels. This information is created by
%                                 hcp_baddata.m pipeline.
%                .badsegmfile: This is the file containing information about bad segments. This information is created by
%                                 hcp_baddata.m pipeline.
%                .icainfofile: This is the file containing information about artifactual Independent Components in the data. This
%                                 information is created hcp_icaclass.m pipeline
%                .badsegmode:  This variable defines if trials containing bad segments will be removed in full or the bad segments will be replaced by NANs.
%                              'remfull' for remove full trial or 'repnan'
%                              for replacing with NANs.  This field is set in the alltrialdefparams_*.m scripts.
%
%                .montage:     This is the montage containing the emg channels. This
%                                variable is constructed by the hcp_exgmontage.m function.
%                                In .labelnew subfield , the expected emg channels names are  expected to be
%                                {'EMG_LH','EMG_LF','EMG_RH','EMG_RF'}
%                .lineFreq:    Numerical array that contain the frequencies
%                                of line current to be filtered out i.e. [60 120].
%                .outputfile:   The filename of the data file where the cleaned data will be saved.
%
%
% OUTPUT
%-----------------------------------------------------------
% cleanTrl : %   This is a numerical matrix similar to the input inCfg.trl trials definition matrix described above.
%                It has 3 columns and number of rows equal to number of CLEAN trials that remained after the cleaning performed.
%                The first column is the start sample, the second is the end sample and time offset of the start of the trial
%                relative to the 0 reference point.
%-------------------------------------------------------------

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.
% Edited by Yasir

datafile = inCfg.datafile;
icainfofile = inCfg.icainfofile;
trl = inCfg.trl;

outputfile = inCfg.outputfile;
badsegmode = inCfg.badsegmode; % tmpCfg.badsegmode = 'remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. when each trial contins an entire block)
montage = inCfg.montage;       % EEG and or EMG (for both one can use E*)
lineFreq = inCfg.lineFreq;

badchanfile = inCfg.badchanfile;
badsegmfile = inCfg.badsegmfile;

useTrlStd = 0;

cfg = [];
cfg.trl = trl;
cfg.datafile = datafile;
cfgDefTr = ft_definetrial(cfg);
cfgDefTr.dataformat = '4d';
cfgDefTr.headerformat = '4d';
cfgDefTr.demean = 'no';
dataRaw = ft_preprocessing(cfgDefTr);

% remove the Supine balancing coefficients, because it's likely incorrect.
% the call to ft_datatype_sens is needed because ft_apply_montage strips
% the chanunit and chantype for some reason
dataRaw.grad = ft_datatype_sens(ft_apply_montage(dataRaw.grad, dataRaw.grad.balance.Supine, 'inverse', 'yes', 'keepunused', 'yes'));

origFs = dataRaw.fsample;
% =====================
% -- Define new sampling frequency to 500 Hz
if origFs>2000
    newFs = origFs/4;
elseif origFs>1000
    newFs = origFs/2;
else
    newFs = origFs;
end
% ========================

elecChans = montage.labelnew';
emgExpChans = {'EMG_LH'; 'EMG_RH';'EMG_LF'; 'EMG_RF'};
[indElec, dummyInd] = match_str(elecChans, emgExpChans);

if ~isempty(indElec)
    tmpmontage = montage;
    tmpmontage.tra = montage.tra(indElec, :);
    hasELEC = 1;
    emgChans = elecChans(indElec);
    tmpmontage.labelnew = emgChans;
    icfg = []; icfg.channel = tmpmontage.labelorg;
    dataRawELECorg = ft_selectdata(icfg, dataRaw);
    dataRawELECnew = ft_apply_montage(dataRawELECorg, tmpmontage);clear dataRawELECorg;
    icfg = []; icfg.channel = emgChans;
    dataRawELEC = ft_selectdata(icfg, dataRawELECnew);clear dataRawELECnew;
    
else
    hasELEC = 0;
    dataRawELEC = [];
end

% ===========================================

hcp_read_matlab(icainfofile, 'comp_class');
hcp_read_ascii(badchanfile);% badchannel
hcp_read_ascii(badsegmfile); % badsegment

badChannels = badchannel.all;
badSegments = badsegment.all;

% ==================================================
inData = dataRaw; clear dataRaw;

% ===============================================
if ~isempty(badChannels)
    badChannelsLabels = cellfun(@(x) ['-', x], badChannels, 'UniformOutput', false);
else
    badChannelsLabels = [];
end

% -------------------------------------------------------------------
selChan = {'MEG', 'MEGREF'};
if ~isempty(badChannelsLabels),
    selChan = [selChan, badChannelsLabels];
end
% -------------------------------------------------------------------
cfg = [];
cfg.channel = selChan;
dataClean0 = ft_selectdata(cfg,inData); clear inData;

% =======================================
cfg = [];
cfg.channel = {'MEG', 'MEGREF'};
dataNEURO = ft_selectdata(cfg, dataClean0);
clear dataClean0;

%% Take out beginning and end
% The beginning and end of every rest scan is marked bad in a way that
% ensures there is always 5 minutes at the middle
% For task, there is alway a part that has to be taken out at start which
% marks before the first trigger
if badSegments(1,1) == 1 % Just to be sure
    beginend = badSegments([1], :); % Removed end that was at hcp_repair
    middle = badSegments(2:(end), :); % end-1 to end
    % Reject
    cfg = [];
    cfg.artfctdef.reject = 'complete';
    cfg.artfctdef.all.artifact = beginend;
else
    middle = badSegments;
end

dataNEURO = ft_rejectartifact(cfg, dataNEURO);
%% Reject bad segments
% =======================================

% =======================================
% -- Reject bad segments
if strcmp(badsegmode, 'remfull')
    cfg = [];
    cfg.artfctdef.reject = 'complete';
    cfg.artfctdef.all.artifact = badSegments;
    
elseif strcmp(badsegmode, 'partial')
    cfg = [];
    cfg.artfctdef.reject = 'partial';
    cfg.artfctdef.all.artifact = badSegments;    
    
elseif strcmp(badsegmode, 'repnan')
    disp('Trial structure preserved. Bad segments replaced with nan');
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.all.artifact = middle;
    
    origDataNEURO = dataNEURO; % rmfield(dataNEURO, 'trial');
%     if hasELEC
%         origDataELEC = dataRawELEC; % rmfield(dataRawELEC, 'trial');
%     end
    
end

if strcmp(badsegmode, 'remfull')|strcmp(badsegmode, 'partial')
    
    dataCleanNEURO1 = ft_rejectartifact(cfg, dataNEURO);
    if hasELEC
        dataCleanELEC1 = ft_rejectartifact(cfg, dataRawELEC);
    end
    
elseif strcmp(badsegmode, 'repnan')
    dataCleanNEURO1 = ft_rejectartifact(cfg, dataNEURO);
%     dataCleanNEURO1 = rmfield(dataNEURO, {'trial', 'time', 'trialinfo', 'sampleinfo'});
%     dataCleanNEURO1.trial = {};
%     dataCleanNEURO1.time = {};
%     dataCleanNEURO1.sampleinfo = [];
%     
%     for iTr = 1:length(dataNEURO.trial)
%         cfg2 = []; cfg2.trials = false(1, length(dataNEURO.trial));
%         cfg2.trials(iTr) = true;
%         tmptrialdata = ft_selectdata(cfg2, dataNEURO);
%         try
%             tmpcleantrialdata = ft_rejectartifact(cfg, tmptrialdata);
%             
%             
%             dataCleanNEURO1.trial = [dataCleanNEURO1.trial tmpcleantrialdata.trial];
%             dataCleanNEURO1.time = [dataCleanNEURO1.time tmpcleantrialdata.time];
%             dataCleanNEURO1.sampleinfo = [dataCleanNEURO1.sampleinfo ; tmpcleantrialdata.sampleinfo];
%         catch
%             disp(['It seems that trial: ', num2str(iTr), ' falls entirely within a artifactual period']);
%         end
%     end
    
    if hasELEC
        % God help us
        dataCleanELEC1 = ft_rejectartifact(cfg, dataRawELEC);
    %     dataCleanELEC1 = rmfield(dataRawELEC, {'trial', 'time', 'trialinfo', 'sampleinfo'});
    %     dataCleanELEC1.trial = {};
    %     dataCleanELEC1.time = {};
    %     dataCleanELEC1.sampleinfo = [];
        
    %     for iTr = 1:length(dataRawELEC.trial)
    %         tmptrialdata = ft_selectdata(dataRawELEC, 'rpt', iTr);
    %         try
    %             tmpcleantrialdata = ft_rejectartifact(cfg, tmptrialdata);
                
    %             dataCleanELEC1.trial = [dataCleanELEC1.trial tmpcleantrialdata.trial];
    %             dataCleanELEC1.time = [dataCleanELEC1.time tmpcleantrialdata.time];
    %             dataCleanELEC1.sampleinfo = [dataCleanELEC1.sampleinfo ; tmpcleantrialdata.sampleinfo];
    %         catch
    %             disp(['It seems tha trial: ', num2str(iTr), ' falls entirely within a artifactual period']);
    %         end
    %     end
    end % if hasELEC
    
end
clear dataRawELEC dataNEURO;
% Save nan trials. From now on, we will continue without them.
save([outputfile, '_includingnans.mat'], 'dataCleanNEURO1', '-v7.3')
if exist('dataCleanELEC1', 'var')
    save([outputfile, '_includingnans.mat'], 'dataCleanELEC1', '-v7.3')
end

% Add nanflags
nanflags = zeros(length(dataCleanNEURO1.trial), 1);
for i = 1:length(dataCleanNEURO1.trial)
    if any(any(isnan(dataCleanNEURO1.trial{i})))
        nanflags(i) = 1;
    end
end

% =======================================
cfg = [];
cfg.demean = 'yes';
cfg.trials = find(~nanflags)'; % No NAN trials
dataCleanNEURO1 = ft_preprocessing(cfg, dataCleanNEURO1);
if hasELEC
    dataCleanELEC1 = ft_preprocessing(cfg, dataCleanELEC1);
end
% if strcmp(badsegmode, 'repnan')
%     origDataNEURO = ft_preprocessing(cfg, origDataNEURO);
%     if hasELEC
%         origDataELEC = ft_preprocessing(cfg, origDataELEC);
%     end
% end
%clear origDataELEC;
% ===================================================
% -- Denoise with Reference Sensors --
%cfg = [];
%denDataNEURO = ft_denoise_pca(cfg, dataCleanNEURO1);
%clear dataCleanNEURO1 ;

denDataNEURO =  dataCleanNEURO1;
clear dataCleanNEURO1 ;
% =======================================
% -- Remove line noise --------------------
for iLine = 1:length(lineFreq)
    cfg = [];
    cfg.detrend = 'no';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = lineFreq(iLine)+[-1 1];
    denDataNEURO = ft_preprocessing(cfg, denDataNEURO);
    if hasELEC
        dataCleanELEC1 = ft_preprocessing(cfg, dataCleanELEC1);
    end
end
dataCleanNEURO2 = denDataNEURO; clear denDataNEURO;
if hasELEC
    dataCleanELEC2 = dataCleanELEC1; clear dataCleanELEC1;
end

% =======================================
%denDataNEURO = dataCleanNEURO2;

% =======================================
% --- identify trials that have high std
% === THIS IS NOT USE BY DEFAULT =======
% ---- Maybe it should be removed
% I removed it --- YÇ
Ntrials = length(dataCleanNEURO2.trial);


denDataNEURO = dataCleanNEURO2; clear dataCleanNEURO2;
% ======================================

% ======================================
badicacomps = comp_class.class.ecg_eog_ic;
if ~isempty(badicacomps)
    cfg = [];
    cfg.component = badicacomps;
    tmpcomp = comp_class;
    tmpcomp.trial{1} = zeros(length(tmpcomp.label), 1); % dummy field
    tmpcomp.time{1} = 0;% dummy field
    dataNEURO = ft_rejectcomponent(cfg, tmpcomp, denDataNEURO); clear denDataNEURO;
else
    dataNEURO = denDataNEURO; clear denDataNEURO;
end

% ===================================================
% -- Denoise with Reference Sensors --
cfg = [];
dataNEURO = ft_denoise_pca(cfg, dataNEURO);



% ==========================================
% Do resampling
% -- Resample to 500 Hz
rscfg = [];
rscfg.detrend = 'no';
rscfg.resamplefs = newFs;
dataNEUROres = ft_resampledata(rscfg, dataNEURO); %clear dataNEURO; % clear data;
if hasELEC
    dataCleanELEC2res = ft_resampledata(rscfg, dataCleanELEC2); clear dataCleanELEC2;% clear data;
end

% if strcmp(badsegmode, 'repnan')
%     
%     onechdatorig=ft_selectdata(origDataNEURO,'channel',1);
%     onechdatorigres = ft_resampledata(rscfg, onechdatorig); clear onechdatorig;
%     
%     resOrigStartSamps=round(origDataNEURO.sampleinfo(:,1)./4);
%     resOrigSampsPerTrl=cellfun(@(x) length(x),onechdatorigres.time);
%     resStartSamps=round(dataNEURO.sampleinfo(:,1)./4);
%     resSampsPerTrl=cellfun(@(x) length(x),dataNEUROres.time);
%     
%     onechdatorigres.sampleinfo=[resOrigStartSamps resOrigStartSamps+resOrigSampsPerTrl'-1];
%     dataNEUROres.sampleinfo=[resStartSamps resStartSamps+resSampsPerTrl'-1];
%     
%     [dataNEUROres, trlNanFlag] = replaceorigwithnans(onechdatorigres, dataNEUROres);
%     if hasELEC
%         dataCleanELEC2res.sampleinfo=[resStartSamps resStartSamps+resSampsPerTrl'-1];
%         dataCleanELEC2res = replaceorigwithnans(onechdatorigres, dataCleanELEC2res); clear dataCleanELEC2res;
%     end
%     clear onechdatorigres;
%     
% end

clear dataNEURO origDataNEURO; % clear data;
% =======================================
% -- Append ELEC and MEG/EEG

if hasELEC
    data = ft_appenddata([], dataNEUROres, dataCleanELEC2res);
    clear dataCleanELEC2res;
else
    data = dataNEUROres;
end
% =======================================
% if strcmp(badsegmode, 'repnan')
%     addTrlColmn = trlNanFlag;
% elseif strcmp(badsegmode, 'remfull')|strcmp(badsegmode, 'partial')
    addTrlColmn = zeros(length(data.trial), 1);
% end
data.trialinfo = [data.trialinfo addTrlColmn];
% =======================================
data.grad = dataNEUROres.grad;
clear dataNEUROres;

% =========================
cleanTrl = data.trialinfo;
varinfo = whos('data');
if (varinfo.bytes/1000000000)>2
    hcp_write_matlab(outputfile, 'data', '-v7.3');
else
    hcp_write_matlab(outputfile, 'data');
end
clear data;

end % function

% ======================================================
% ======================================================
% ======================================================
% ======================================================

function[newdata, trlNanFlag] = replaceorigwithnans(origdata, cleandata)
% The original data is the one before the partial artifact rejection
% The clean data is after the partial artifact rejection
% The function outputs a newdata with the same trial structure as the
% origdata but with the cleandata data in the clean periods and nans
% in the noisyperiosds; It also outputs an array trlNanFlag with length
% equal to the trials in the newdata and flags with 1 the trials that have
% nans and 0's otherwise.
% The original data is the one before the partial artifact rejection
% The clean data is after the partial artifact rejection

newdata = origdata;
origSampinfo = origdata.sampleinfo;
cleanSampinfo = cleandata.sampleinfo;

NtrlOrig = size(origSampinfo, 1);
NtrlClean = size(cleanSampinfo, 1);

Nchans = length(cleandata.label);

trlNanFlag = zeros(NtrlOrig, 1);
for iOrig = 1:NtrlOrig,
    % iOrig
    origStart = origSampinfo(iOrig, 1);
    origEnd = origSampinfo(iOrig, 2);
    tmptrldata = nan(Nchans, origEnd-origStart+1);
    
    for iClean = 1:NtrlClean
        cleanStart = cleanSampinfo(iClean, 1);
        cleanEnd = cleanSampinfo(iClean, 2);
        if (cleanStart == origStart) && (cleanEnd == origEnd)
            tmptrldata = cleandata.trial{iClean};
        else
            [comSamps, indxOrig, indxClean] = intersect([origStart:origEnd], [cleanStart:cleanEnd]);
            if ~isempty(comSamps)
                tmptrldata(:, indxOrig) = cleandata.trial{iClean}(:, indxClean);
            end
        end
        
    end
    if any(isnan(tmptrldata(1, :)))
        trlNanFlag(iOrig) = 1;
    end
    newdata.trial{iOrig} = tmptrldata;
end
newdata.label = cleandata.label;

end % function
