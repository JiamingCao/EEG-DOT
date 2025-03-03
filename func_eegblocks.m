function func_eegblocks(filename)

% Preprocessing: filter, detrend, and reref
cfg = [];
cfg.dataset = filename;
if strcmp(filename, 'Subject104/leftright.bdf') || ~isempty(regexp(filename,'Subject111/*'))
    cfg.channel = [1:64, 257, 258]; % In these files, all 256 channels were recorded!
else
    cfg.channel = [1:66];
end
if ~isempty(regexp(filename,'Subject110/*'))
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [59, 61];
end
cfg.bpfilter = 'yes';
cfg.bpfreq = [1, 50];
cfg.bpfiltord = 4;  % default: butterworth
cfg.polyremoval = 'yes';    % default: 2nd order
cfg.reref = 'yes';
cfg.refchannel = [65, 66];    % avg of mastoid channels   
data = ft_preprocessing(cfg);
newlabel = {'Fp1', 'AF7', 'AF3', 'F1', 'F3', 'F5', 'F7', 'FT7', 'FC5', 'FC3', 'FC1', 'C1', 'C3', 'C5', 'T7', 'TP7', ...
    'CP5', 'CP3', 'CP1', 'P1', 'P3', 'P5', 'P7', 'P9', 'PO7', 'PO3', 'O1', 'Iz', 'Oz', 'POz', 'Pz', 'CPz', ...
    'Fpz', 'Fp2', 'AF6', 'AF4', 'AFz', 'Fz', 'F2', 'F4', 'F6', 'F8', 'FT8', 'FC6', 'FC4', 'FC2', 'FCz', 'Cz', ...
    'C2', 'C4', 'C6', 'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO8', 'PO4', 'O2'}';
data.label(1:64) = newlabel;
% data.hdr.label = newlabel;
fs_original = data.fsample;
% resample
fs_new = 256;   % Hz
cfg = [];
cfg.resamplefs = fs_new;
data = ft_resampledata(cfg, data);
% ICA
cfg = [];
cfg.channel = 1:64;
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
% cfg.method = 'fastica';
comp = ft_componentanalysis(cfg, data);
% reject EOG
eog = mean(data.trial{1}([1, 33, 34], :));
R = corr([eog', comp.trial{1}']);
comp_reject = find(R(1, 2:end) > 0.4);
cfg = [];
cfg.component = comp_reject;
data = ft_rejectcomponent(cfg, comp);
% interpolate bad channels
cfg=[];
cfg.method='template';
cfg.template='biosemi64_neighb.mat';
neighbors=ft_prepare_neighbours(cfg, data);
channel_std = std(data.trial{1}, [], 2);
bad_channels = data.label(channel_std > 50);   % manually determined threshold for "too much fluctuation"

cfg = [];
cfg.method = 'average';
cfg.badchannel = bad_channels;
cfg.neighbours = neighbors;
cfg.layout='biosemi64.lay';
data = ft_channelrepair(cfg, data);
% get triggers
cfg = [];
cfg.dataset = filename;
cfg.trialdef.eventtype = 'STATUS';
cfg.trialdef.eventvalue = {65281, 65282};
events = ft_definetrial(cfg);

trigger1 = round(events.trl(events.trl(:, 4) == 65281, 1) * fs_new/fs_original);
trigger2 = round(events.trl(events.trl(:, 4) == 65282, 1) * fs_new/fs_original);
trigger1(find(diff(trigger1) < 0.4*fs_new)) = [];
trigger2(find(diff(trigger2) < 0.4*fs_new)) = [];
if strcmp(filename, 'Subject101/topbottom.bdf')
    trigger1 = trigger1(end-300+1:end); % fixing an error in subject 101, topbottom
    trigger2 = trigger2(end-300+1:end);
end

%% Block average
EEG = data.trial{1};
if ~isempty(regexp(filename,'Subject116/*'))
    EEG = PCAFilter(EEG', 7)';
else
    EEG = PCAFilter(EEG', 1)';
end
prestim = 0.05; % sec
poststim = 0.75; % sec
prestim_pt = round(prestim*fs_new);
poststim_pt = round(poststim*fs_new);
t_block = linspace(-prestim, poststim, prestim_pt+poststim_pt)*1000; % ms
if ~isempty(regexp(filename,'Subject101/*')) || ~isempty(regexp(filename,'Subject102/*'))
    trigger1 = trigger2;
    trigger1(20:20:end) = [];
    fprintf('Tiggers flipped\n')
end
block = zeros(64, length(trigger1), prestim_pt+poststim_pt);
for ch = 1:64
    for i = 1:length(trigger1)
        tmp = EEG(ch, trigger1(i)-prestim_pt:trigger1(i)+poststim_pt - 1);
        block(ch, i, :) = tmp - mean(tmp(1:prestim_pt));
    end
end

for i=1:length(trigger1)
    if max(max(abs(squeeze(block(:,i,:))))) > 100
        block(:,i,:) = nan;
    end
end
fprintf('Block: %d trials rejected\n', sum(isnan(block(:)))/(size(block, 1)*size(block, 3)));
avg_blk = squeeze(nanmean(block, 2));
% figure,plot(t_block, avg_blk(29,:))

%% save data
tmp = split(filename, '.');
save([tmp{1}, '_processed'],'EEG','avg_blk','bad_channels','block','ch','fs_new','poststim','poststim_pt','prestim', 'prestim_pt','t_block','trigger1','trigger2')
fprintf(['Done processing ', filename, '\n'])
