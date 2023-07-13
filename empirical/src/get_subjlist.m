%% Get subject list that has resting state
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

cd(datafolder)
task = string('Restin');

subjlist = dir;
subjlist = {subjlist.name};
subjs = [];
for i = 3:length(subjlist) % Ignore . and ..
    paths = genpath(subjlist{i});
    if contains(paths, task)
        subjs = [subjs; subjlist{i}];
    end
end

save(char(string(projectfolder) + string('/data/') + task + string('_subjs.mat')), 'subjs') % Sorry...
%% Working Memory
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

cd(datafolder)
task = string('Wrkmem');

subjlist = dir;
subjlist = {subjlist.name};
subjs = [];
for i = 3:length(subjlist) % Ignore . and ..
    paths = genpath(subjlist{i});
    if contains(paths, task)
        subjs = [subjs; subjlist{i}];
    end
end

save(char(string(projectfolder) + string('/data/') + task + string('_subjs.mat')), 'subjs')
disp('EVERYTHING FINISHED!!!')
%% Story-Math
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

cd(datafolder)
task = string('StoryM');

subjlist = dir;
subjlist = {subjlist.name};
subjs = [];
for i = 3:length(subjlist) % Ignore . and ..
    paths = genpath(subjlist{i});
    if contains(paths, task)
        subjs = [subjs; subjlist{i}];
    end
end

save(char(string(projectfolder) + string('/data/') + task + string('_subjs.mat')), 'subjs')
disp('EVERYTHING FINISHED!!!')
%% MotorTask
cd /BICNAS2/ycatal/sourcereconstruction
set_paths

cd(datafolder)
task = string('Motort');

subjlist = dir;
subjlist = {subjlist.name};
subjs = [];
for i = 3:length(subjlist) % Ignore . and ..
    paths = genpath(subjlist{i});
    if contains(paths, task)
        subjs = [subjs; subjlist{i}];
    end
end

save(char(string(projectfolder) + string('/data/') + task + string('_subjs.mat')), 'subjs')
disp('EVERYTHING FINISHED!!!')