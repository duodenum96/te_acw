%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by duodenum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup the execution environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl software;

% ensure that the time and date of execution are not stored in the provenance information
% global ft_default
% ft_default.trackcallinfo = 'no';
% 
% % allow the user to specify the path where additional data is present, e.g. the channel layout or anatomy files
if exist('path', 'var')
    addpath(path)
end

if ~exist('filename', 'var')
    error('filename should be specified')
end
% 
% % the filename is assumed to be something like
% % 'rawdatadir/Phase1MEG/Subjects/CP10018/Experiments/CP10018_MEG/Scans/1-Rnoise_MNN_V1/Resources/4D/c,rfDC'
% tok = tokenize(filename, '/');
% 
% if ~exist('subjectid', 'var')
%     subjectid = tok{end-7};
% end
% 
% if ~exist('experimentid', 'var')
%     experimentid = tok{end-5};
% end
% 
% if ~exist('scanid', 'var')
%     scanid = tok{end-3};
% end
% 
% if ~exist('pipelinedatadir', 'var')
%     pipelinedatadir = hcp_pathdef;
% end


resultprefix = sprintf('%s_%s_%s-%s', subjectid, experimentid, scanid, task);
% resultprefix = sprintf('%s_%s', experimentid, scanid);

% change to the location of the processed data (input and output)
cd(pipelinedatadir)

% hcp_check_pipelineoutput('baddata', subjectid, scanid, task);
% hcp_check_pipelineoutput('icaclass', subjectid, scanid, task);
% hcp_check_pipelineoutput('icaclass_qc', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
% hcp_check_pipelineoutput('anatomy', subjectid, scanid, task);


ver('megconnectome')

% print the value of all local variables to screen for provenance
% w = whos;
% w = {w.name};
% w = setdiff(w, {'w', 'ans'});
% for i=1:length(w)
%     fprintf(hcp_printstruct(w{i}, eval(w{i})));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hcp_read_matlab([pipelinedatadir, '/Restin/icaclass/', resultprefix, '_icaclass_vs.mat'])
hcp_read_matlab([pipelinedatadir, '/', task, '/icaclass/', resultprefix '_icaclass_vs.mat'])

% read the source and volume conduction model from current dir with
% outputs of previous pipelines
% No.

hcp_read_matlab([pipelinedatadir, '/anatomy/', subjectid, '_MEG_anatomy_sourcemodel_2d']);
sourcemodel2d=ft_convert_units(sourcemodel2d, 'cm');
sourcemodel2d.inside = 1:size(sourcemodel2d.pos,1);
sourcemodel2d.outside = [];
sourcemodelsubj = sourcemodel2d;

hcp_read_matlab([pipelinedatadir, '/anatomy/', subjectid, '_MEG_anatomy_headmodel']);
headmodel = ft_convert_units(headmodel, 'cm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grad = ft_read_sens(filename);
gradBalanced = grad;
gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
grad=gradBalanced;
grad = ft_convert_units(grad, 'cm');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the component data in order for ft_sourceanalysis to be able to
% swallow it
mixing   = comp_class.topo;
channels = comp_class.topolabel;
% normalisation of the topographies
for i = 1:size(mixing, 2)
  val(i) = 0.01*max(abs(mixing(:, i)));
  mixing(:, i) = mixing(:, i)/val(i);
end

% create a 'timelock' structure
tlck = [];
tlck.label = channels;
tlck.cov = eye(numel(tlck.label)); 
tlck.time=1;
tlck.grad = grad;
tlck.dimord = 'chan_time';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the forward solution

cfg = [];
cfg.headmodel = headmodel;
cfg.sourcemodel = sourcemodelsubj;
cfg.grad = grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;
gridLF = ft_prepare_leadfield(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do an MNE with different regularisation for each component

% this parameter is hard-coded here
noise_level = 8;

% specify the static part of the cfg for the source reconstruction
cfg               = [];
cfg.method        = 'mne';
cfg.grid          = gridLF;
cfg.vol           = headmodel;
cfg.channel       = channels;
cfg.mne.prewhiten = 'yes';
cfg.mne.noisecov  = eye(numel(channels))*noise_level;
cfg.mne.keepfilter = 'yes';
cfg.mne.scalesourcecov = 'yes';

% loop over components, due to component-specific regularisation
for i=1:size(mixing,2)

  % use the channel-level topography of the current component
  tlck.avg = mixing(:,i);
  
  % estimate the snr of the current component
  cfg.mne.snr = sqrt(mean((mixing(:,i)-mean(mixing(:,i))).^2))/noise_level;
  noisevec(i) = cfg.mne.snr;
  
  tmp = ft_sourceanalysis(cfg, tlck);
  inside_indices = find(tmp.inside(:))';
  if i==1
    % create the output source structure in the first iteration
    source=tmp;
  else
    % concatenate the reconstructed source level topography to the previously computed ones
    for k = 1:numel(inside_indices)
      source.avg.mom{inside_indices(k)} = cat(2,source.avg.mom{inside_indices(k)}, tmp.avg.mom{inside_indices(k)});
    end
    source.avg.pow = horzcat(source.avg.pow,tmp.avg.pow);
  end
end

% add some relevant fields
source.val  = val;
source.time = 1:size(mixing,2);
source.snr  = noisevec;
if isfield(sourcemodelsubj,'tri'),            source.tri            = sourcemodelsubj.tri;            end
if isfield(sourcemodelsubj,'brainstructure'), source.brainstructure = sourcemodelsubj.brainstructure; end
if isfield(sourcemodelsubj,'brainstructurelabel'), source.brainstructurelabel = sourcemodelsubj.brainstructurelabel; end

% save the data as a mat-file
outfolder = ['/BICNAS2/ycatal/sourceparcellated/', subjectid, '/', task, '/icamne/'];
if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end
save([outfolder, resultprefix,'_icamne'], 'source')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting and QC is done here

% this is for plotting
% cfgtopo =[];
% cfgtopo.grad=grad;
% cfgtopo.zlim='maxabs';
% cfgtopo.component=[];
% cfgtopo.parameter = 'topo';
% cfgtopo.comment = 'no';
% cfgtopo.colorbar = 'yes';
% cfgtopo.layout='4D248.mat';
% tmpclass=comp_class;
% tmpclass.trial{1}(1,1)=0;
% tmpclass.time{1}(1,1)=0;
% tmpclass.topo=mixing;
% 
% X = 29.7;                  %# A3 paper size
% Y = 14.35;                  %# A3 paper size
% xMargin = 1;               %# left/right margins from page borders
% yMargin = 1;               %# bottom/top margins from page borders
% xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
% ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
% 
% 
% cfg = [];
% cfg.funparameter = 'avg.pow'; 
% cfg.method = 'surface';
% cfg.funcolormap='jet';
% cfg.interactive = 'no';
% 
% tmp = source;
% for k = 1:numel(source.time)
%     tmp.avg.pow = source.avg.pow(:, k);
%     tmp.time=source.time(k);
%     maxabs=max(max(max(tmp.avg.pow)));
%     imgname = [resultprefix '_icamne_' num2str(k) '.png'];
%     
%     ft_plot_mesh(tmp,'vertexcolor',tmp.avg.pow)
%     view(-90,90)
%     
%     h1=gcf;
%     set(gca, 'XTickLabel',[], 'YTickLabel',[], ...
%         'Units','normalized', 'Position',[0 0 1 1])
%     
%     %# figure size on screen (50% scaled, but same aspect ratio)
%     set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize])
%     
%     %# figure size printed on paper
%     set(gcf, 'visible', 'on')
%     set(gcf, 'PaperUnits','centimeters')
%     set(gcf, 'PaperSize',[X Y])
%     set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
%     set(gcf, 'PaperOrientation','portrait')
%     
%     fax(1)=gca;
%     set(fax(1),'position',[0.25 0.1 0.95 0.8]);
%     fax(2)= axes('position',[0 0.1 0.5 0.8]);
%     
%     cfgtopo.component=k;
%     ft_topoplotIC(cfgtopo, tmpclass);
%     
%     hcp_write_figure(imgname, h1)
%     close(h1)
% end
% 
% % ensure that the expected output files were created
% hcp_check_pipelineoutput('icamne', 'subject', subjectid, 'experiment', experimentid, 'scan', scanid);
