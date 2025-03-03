clear;
addpath('~/Documents/MATLAB/fieldtrip');
ft_defaults;

% load('groupeeg_topbottom.mat')
% load('groupnirs_topbottom.mat')
load('groupeeg_leftright.mat')
load('groupnirs_leftright.mat')
% load('groupnirs_leftright_project.mat')
load('templatebrain/eegmodel8196.mat')
% load('templatebrain/eegmodel20484.mat')

%% EEG-NIRS
% nirspow=abs(recon_groupbeta);
% nirspow = abs(Beta_Hb_proj);
nirspow = abs(avg_recon1);
thresh = 0;
nirspow_norm = nirspow.*(nirspow > thresh) / max(nirspow);
Qn_eeg = {speye(size(L, 1))};
% Qp_eeg = {0.1*speye(size(L,2))+0.9*spdiags(double(nirspow_norm>0.1),0,size(L, 2), size(L, 2))};
Qp_eeg = {spdiags(1-exp(-(nirspow_norm + 0.1)/1), 0, size(L, 2), size(L, 2))};
% Qp_eeg = {spdiags(1-exp(-(nirspow_norm + 0.1)/0.5), 0, size(L, 2), size(L, 2))};

% [~, eegrecon21, inv_op21] = REML(100*avg1(:, 58), L, recon1, Qn_eeg, Qp_eeg, 500);  % Estimate inverse operator using N170 peak
% [~, eegrecon22, inv_op22] = REML(100*avg2(:, 58), L, recon2, Qn_eeg, Qp_eeg, 500);
% 
% alleeg21 = inv_op21*avg1;
% alleeg22 = inv_op22*avg2;
blkeeg2 = zeros(size(L,2), length(t_block));
for i=1:length(t_block)
    [~, blkeeg2(:,i)] = REML(100*blk_avg(:, i), L, [], Qn_eeg, Qp_eeg, 500);
end

% Twomey
% C=diag(var(blk_avg(:,1:25),0,2));
% lambda = 1e-2;
% tmp1=inv(L'*inv(C)*L+lambda*eye(8196));
% tmp2=L'*inv(C)*blk_avg;
% blkeeg3=tmp1*(tmp2+lambda*blkeeg2/100);

% figure,ft_plot_mesh(source, 'vertexcolor', eegrecon21, 'colormap', redblue);view(0,0);title('eeg nirs 1')
% caxis([-max(abs(eegrecon21)), max(abs(eegrecon21))])
% figure,ft_plot_mesh(source, 'vertexcolor', eegrecon22, 'colormap', redblue);view(0,0);title('eeg nirs 2')
% caxis([-max(abs(eegrecon22)), max(abs(eegrecon22))])

%%
% figure;
% set(gcf, 'Position', [150,160,1600,800])
% % v = VideoWriter('eeg_topbottom_prior.avi');
% v = VideoWriter('eeg_leftright_prior.avi');
% v.FrameRate = 12;
% open(v);
% for i=1:length(t_epoch)
%     subplot(1,2,1), cla;
%     ft_plot_mesh(source, 'vertexcolor', alleeg21(:,i), 'colormap', redblue);view(0,0);
%     camlight
%     caxis([-325, 325]);
%     title(['t=', num2str(t_epoch(i)), ' ms']);
%     subplot(1,2,2), cla;
%     ft_plot_mesh(source, 'vertexcolor', alleeg22(:,i), 'colormap', redblue);view(0,0);
%     camlight
%     caxis([-290, 290])
%     title(['t=', num2str(t_epoch(i)), ' ms']);
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% %     pause(0.05);
% end
% close(v)

%%
figure;
set(gcf, 'Position', [150,160,800,800])
% v = VideoWriter('eeg_topbottom_blk_prior.avi');
v = VideoWriter('eeg_leftright_blk_prior.avi');
v.FrameRate = 12;
open(v);
for i=1:length(t_block)
    cla;
    ft_plot_mesh(source, 'vertexcolor', avg_prior(:,i), 'colormap', redblue);view(0,0);
    camlight headlight
    caxis([-8e4, 8e4]);
    title(['t=', num2str(t_block(i)), ' ms']);
    frame = getframe(gcf);
    writeVideo(v, frame);
%     pause(0.05);
end
close(v)

%%
figure;
set(gcf, 'Position', [150,160,800,800])
% v = VideoWriter('eeg_topbottom_blk_prior.avi');
v = VideoWriter('eeg_leftright_blk_prior_soft.avi');
v.FrameRate = 12;
open(v);
for i=1:length(t_block)
    cla;
    ft_plot_mesh(source, 'vertexcolor', blkeeg3(:,i), 'colormap', redblue);view(0,0);
    camlight
    caxis([-300, 300]);
    title(['t=', num2str(t_block(i)), ' ms']);
    frame = getframe(gcf);
    writeVideo(v, frame);
%     pause(0.05);
end
close(v)
