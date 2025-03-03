clear;
addpath('~/Documents/MATLAB/fieldtrip');
ft_defaults;

for subject = 110:119
    for experiment = {'leftright', 'topbottom'}
        filename = ['Subject', num2str(subject), '/', experiment{1}, '.bdf'];
        if ~exist(filename, 'file')
            fprintf([filename, ' doesn''t exist. Skipping\n'])
            continue
        elseif strcmp(filename, 'Subject103/topbottom.bdf')
            continue
        else
            fprintf([filename, '\n']);
            func_eegblocks(filename);
        end
    end
end
