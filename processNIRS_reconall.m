addpath(genpath('~/Documents/MATLAB/easyh5'))
addpath(genpath('~/Documents/MATLAB/jsnirfy'))
addpath('~/Documents/MATLAB/fieldtrip');
addpath('fromNIRSToolbox')
ft_defaults;
clear;

source = ft_read_headshape('cortex_8196.surf.gii');
load('templatebrain/nirsmodel8196_full.mat')
% load('fsaverage5/sourcemodel.mat')
% load('fsaverage5/nirsmodel_fs_full.mat')

for subject = 101:119
    for experiment = {'leftright','topbottom'}
        folder = ['Subject', num2str(subject), '/', experiment{1}];
        if ~exist(folder, 'dir')
            fprintf([folder, ' doesn''t exist. Skipping\n'])
            continue
        elseif strcmp(folder, 'Subject112/topbottom')
            continue
        else
            fprintf([folder, '\n']);
            func_nirsrecon(folder, source, jacobian_full, nirs_mesh);
        end
    end
end
