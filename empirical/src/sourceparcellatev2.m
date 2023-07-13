%% Source level parcellated results
% Inspired by Soren's code
%% Import files

hcp_read_matlab([pipelinedatadir, '/anatomy/', subjectid, '_MEG_anatomy_sourcemodel_2d']); % sourcemodel2d
sourcemodel2d = ft_convert_coordsys(sourcemodel2d, 'mni');

% Load stuff
tempfolder = '/HCP/MEG/templates';
atlas = ft_read_atlas([tempfolder, '/atlas_MMP1.0_4k.mat']);
atlas = ft_convert_units(atlas,'cm');

if string(task) == string('Restin')
    datafolder = ['/BICNAS2/ycatal/sourceparcellated/', subjectid, '/', task, '/rmegpreproc'];
    data = load([datafolder, '/', subjectid, '_MEG_', scanid, '-', task, '_rmegpreproc.mat']); data = data.data;
else
    datafolder = ['/BICNAS2/ycatal/sourceparcellated/', subjectid, '/', task, '/tmegpreproc'];
    data = load([datafolder, '/', subjectid, '_MEG_', scanid, '-', task, '_tmegpreproc_', taskextension, '.mat']); data = data.data;
end

load(['/BICNAS2/ycatal/sourceparcellated/', subjectid, '/', task, '/icamne/', subjectid, '_MEG_', scanid, '-', task, '_icamne.mat'])
% load([pipelinedatadir, '/', task, '/icamne/', subjectid, '_MEG_', scanid, '-', task, '_icamne.mat']) % source

% Select only the relevant channels
cfg = []; cfg.channel = {'MEG'};
data = ft_selectdata(cfg,data);

source.avg.filter2 = [];
for i = 1:length(data.label)
    chanindx = find(strcmp(data.label(i), source.avg.label));
    for j = 1:length(source.avg.filter)
        source.avg.filter2{j}(:, i) = source.avg.filter{j}(:, chanindx);
    end
end
source.avg.filter = []; source.avg.filter = source.avg.filter2; source.avg.filter2 = [];

% If continuous, epoch into arbitrary 2-second segments
if length(data.trial) == 1
    cfg = []; cfg.event = 1:2*data.fsample:length(data.time{1}); cfg.epoch = [0 (2*data.fsample)-1];
    cfg.event(end) = [];
    data = ft_epoch(cfg,data);
end
atlas.pos = source.pos; % Visually checked just to be sure, they're in the same space / coordinates
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'parcellation';
interpsource = ft_sourceinterpolate(cfg, atlas, source);
interpsource.avg = source.avg;

cfg = [];
cfg.parcellation = 'parcellation';
roidata = ft_virtualchannel(cfg, data, interpsource, atlas);

outfoldersrc = [outfolder '/source/'];
if ~exist(outfoldersrc, 'dir')
    mkdir(outfoldersrc)
end
if exist('taskextension', 'var')
    resultprefix = sprintf('%s_%s_%s-%s_%s', subjectid, experimentid, scanid, task, taskextension);
else
    resultprefix = sprintf('%s_%s_%s-%s', subjectid, experimentid, scanid, task);
end

save([outfoldersrc, resultprefix,'_glasser'], 'roidata')
