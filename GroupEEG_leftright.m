clear;
load('templatebrain/eegmodel8196.mat')
% load('fsaverage5/eegmodel_fs.mat')
Subjects = dir('Subject*');
% epochs1 = [];
% epochs2 = [];
block = [];
bad_nirs = [4,9,15];
for s = 1:length(Subjects)
    if ~exist([Subjects(s).name, '/leftright_processed.mat'], 'file')
        continue
    end
    if any(s == bad_nirs)
        continue
    end
    tmp = load([Subjects(s).name, '/leftright_processed']);
    block = cat(2, block, tmp.block);
end

% t_epoch = tmp.t_epoch;
t_block = tmp.t_block;
elecpos = headmodel.elec.elecpos;

% avg1 = squeeze(nanmean(epochs1, 2));
% avg2 = squeeze(nanmean(epochs2, 2));
blk_avg = squeeze(nanmean(block, 2));

alpha = 5e-7;
[~,invop]=tikhonov(L, alpha);
% recon1=invop*avg1(:,58);    % N170 for demonstration
% recon2=invop*avg2(:,58);
% figure,ft_plot_mesh(source, 'vertexcolor', recon1, 'colormap', redblue);view(0,0);
% caxis([-max(abs(recon1)), max(abs(recon1))]), title('left wedge')
% figure,ft_plot_mesh(source, 'vertexcolor', recon2, 'colormap', redblue);view(0,0);
% caxis([-max(abs(recon2)), max(abs(recon2))]), title('right wedge')
% 
% alleeg1 = invop*avg1;
% alleeg2 = invop*avg2;
blkeeg = invop*blk_avg;
save('groupeeg_leftright.mat', 'alpha', 'blk_avg', 'blkeeg', 'block', 'elecpos', 'invop', 't_block');

%%
% figure;
% set(gcf, 'Position', [150,160,1600,800])
% v = VideoWriter('eeg_leftright.avi');
% v.FrameRate = 12;
% open(v);
% for i=1:length(t_epoch)
%     subplot(1,2,1), cla;
%     ft_plot_mesh(source, 'vertexcolor', alleeg1(:,i), 'colormap', redblue);view(0,0);
%     camlight
%     caxis([-740, 740]);
%     title(['t=', num2str(t_epoch(i)), ' ms']);
%     subplot(1,2,2), cla;
%     ft_plot_mesh(source, 'vertexcolor', alleeg2(:,i), 'colormap', redblue);view(0,0);
%     camlight
%     caxis([-510, 510])
%     title(['t=', num2str(t_epoch(i)), ' ms']);
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% %     pause(0.05);
% end
% close(v)

%%
figure;
set(gcf, 'Position', [150,160,800,800])
v = VideoWriter('eeg_leftright_blk.avi');
v.FrameRate = 12;
open(v);
for i=1:length(t_block)
    cla;
    ft_plot_mesh(source, 'vertexcolor', avg_noprior(:,i), 'colormap', redblue);view(0,0);
    camlight headlight
    caxis([-200, 200]);
    title(['t=', num2str(t_block(i)), ' ms']);
    frame = getframe(gcf);
    writeVideo(v, frame);
%     pause(0.05);
end
close(v)
