%% Set paths for fieldtrip, megconnectome etc.
clear, clc, close all
% restoredefaultpath
ftpath = '/BICNAS2/ycatal/fieldtrip_20221210';
addpath(ftpath)
ft_defaults
projectfolder = '/BICNAS2/ycatal/sourcereconstruction';
cd(projectfolder)
addpath(genpath(projectfolder))
addpath(genpath('/BICNAS2/ycatal/megconnectome-3.0'))

datafolder = '/HCP/MEG/megdata';
