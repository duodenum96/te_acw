function hcp_repair(inCfg)
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
% Modified by Yasir
% Original function: hcp_extract_allfromrun
% Modified to repair instead of remove / replace w/ nans

datafile = inCfg.datafile;
icainfofile = inCfg.icainfofile;
trl = inCfg.trl;

outputfile = inCfg.outputfile;
badsegmod = inCfg.badsegmode; % tmpCfg.badsegmode = 'remfull'; % This is used by tmeg_preproc to know if trials should be fully removed, or bad segments should be replaced by nans(i.e. when each trial contins an entire block)
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
    dataRawELECorg = ft_selectdata(dataRaw, 'channel', tmpmontage.labelorg);
    dataRawELECnew = ft_apply_montage(dataRawELECorg, tmpmontage);clear dataRawELECorg;
    dataRawELEC = ft_selectdata(dataRawELECnew, 'channel', emgChans);clear dataRawELECnew;
    
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
dataClean0 = ft_selectdata(cfg, inData); clear inData
% dataClean0 = ft_selectdata(inData, 'channel', selChan); clear inData;

% =======================================
cfg = [];
cfg.channel = {'MEG', 'MEGREF'};
dataNEURO = ft_selectdata(cfg, dataClean0);
% dataNEURO = ft_selectdata(dataClean0, 'channel', {'MEG', 'MEGREF'});
clear dataClean0;
% =======================================

% =======================================
%% Take out beginning and end
% The beginning and end of every rest scan is marked bad in a way that
% ensures there is always 5 minutes at the middle
% Define beginning and end and also middle for next step
beginend = badSegments([1, end], :);
middle = badSegments(2:(end-1), :);
% Reject
cfg = [];
cfg.artfctdef.reject = 'complete';
cfg.artfctdef.all.artifact = beginend;

dataCleanMiddle = ft_rejectartifact(cfg, dataNEURO);
clear dataNEURO

% Repair
if ~isempty(middle)
    disp('Bad segments detected, will interpolate')

    % Put NaNs in bad segments (rejection)
    cfg = [];
    cfg.artfctdef.reject = 'nan';
    cfg.artfctdef.all.artifact = middle;
    dataClean = ft_rejectartifact(cfg, dataCleanMiddle); clear dataCleanMiddle
    % Put NaNFlag for rejected trials
    nanflags = zeros(length(dataClean.trial), 1);
    for i = 1:length(dataClean.trial)
        if any(any(isnan(dataClean.trial{i})))
            nanflags(i) = 1;
        end
    end
    % Interpolate NaNs
    catdata = cat(2, dataClean.trial{:});
    catdata = fillmissing(catdata', 'linear')';
    
    % Put the catenated data back into ft struct
    trllength = size(dataClean.trial{1}, 2);
    dataCleanInterp = dataClean;
    for i = 1:length(dataClean.trial)
        dataCleanInterp.trial{i} = catdata(:, ((i-1)*trllength+1):(i*trllength));
    end
    clear dataClean catdata
else
    nanflags = zeros(length(dataCleanMiddle.trial), 1);
    dataCleanInterp = dataCleanMiddle; clear dataCleanMiddle % Just to avoid any name mismatch
end
%% Other preprocessing steps
%% Demean
cfg = [];
cfg.demean = 'yes';
dataCleanNEURO1 = ft_preprocessing(cfg, dataCleanInterp); 
clear dataCleanInterp
%% Remove Line Noise
for iLine = 1:length(lineFreq)
    cfg = [];
    cfg.detrend = 'no';
    cfg.bsfilter = 'yes';
    cfg.bsfreq = lineFreq(iLine)+[-1 1];
    dataCleanNEURO2 = ft_preprocessing(cfg, dataCleanNEURO1);
end
clear dataCleanNEURO1
%% Remove bad ICA components
badicacomps = comp_class.class.ecg_eog_ic;
if ~isempty(badicacomps)
    cfg = [];
    cfg.component = badicacomps;
    tmpcomp = comp_class;
%     tmpcomp.trial{1} = zeros(length(tmpcomp.label), 1); % dummy field
%     tmpcomp.time{1} = 0;% dummy field
    dataNEURO = ft_rejectcomponent(cfg, tmpcomp, dataCleanNEURO2); 
else
    dataNEURO = dataCleanNEURO2; 
end
clear dataCleanNEURO2
%% Denoise w/ reference sensors
cfg = [];
dataNEURO2 = ft_denoise_pca(cfg, dataNEURO);
clear dataNEURO
%% Resample to 500 Hz
rscfg = [];
rscfg.detrend = 'no';
rscfg.resamplefs = newFs;
data = ft_resampledata(rscfg, dataNEURO2); 
clear dataNEURO2
%% Add nanflags and save
data.nanflag = nanflags;
varinfo = whos('data');
if (varinfo.bytes/1000000000)>2,
    hcp_write_matlab(outputfile, 'data', '-v7.3');
else
    hcp_write_matlab(outputfile, 'data');
end
end % function
