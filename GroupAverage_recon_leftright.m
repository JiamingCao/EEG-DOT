addpath('~/Documents/MATLAB/fieldtrip');
ft_defaults;

clear;
% source = ft_read_headshape('cortex_8196.surf.gii');
% load('templatebrain/nirsmodel8196_full.mat')
load('fsaverage5/sourcemodel.mat')
load('fsaverage5/nirsmodel_fs_full.mat')
Subjects = dir('Subject*');
%%
dOD_blk = [];
allbeta = [];
allrecon1 = [];
allrecon2 = [];
bad_subject=[];
for s = 1:19
    try
        tmp = load([Subjects(s).name, '/leftright/processed_recon4_fs']);
    catch
        continue
    end
%     bad = load([Subjects(s).name, '/leftright/bad']);
%     dOD_blk = cat(3, dOD_blk, tmp.dOD_blk);
    nchannels = tmp.nchannels;
    beta = nan(2*nchannels,1);
    bad = tmp.bad;
    if length(bad)>30
        bad_subject = [bad_subject, s];
    end
    good = setdiff(1:2*nchannels, [bad;bad+nchannels]);
    beta(good) = tmp.beta;
    allbeta = [allbeta, beta];
    allrecon1 = [allrecon1, tmp.recon1];
    allrecon2 = [allrecon2, tmp.recon2];
end

% bad_subject=[4,9,13];
% bad_subject = [2,6];
% bad_subject=[2];
allrecon1(:,bad_subject) = [];
allrecon2(:,bad_subject) = [];
% blk_avg = mean(dOD_blk, 3, 'omitnan');
allbeta(:,bad_subject) = nan;
avg_beta = nanmean(allbeta, 2);
% avg_recon1 = nanmedian(allrecon1, 2);
% avg_recon2 = nanmedian(allrecon2, 2);
avg_recon1 = nanmean(allrecon1, 2);
avg_recon2 = nanmean(allrecon2, 2);

%% Plotting
figure,ft_plot_mesh(source, 'vertexcolor', avg_recon1/max(abs(avg_recon1)), 'colormap', redblue);view(0,0);title('\DeltaHbO')
caxis([-0.8,0.8]); camlight headlight
saveas(gcf, 'Figures/HbO.fig')
print Figures/HbO.png -dpng -r300
hbo_masked = avg_recon1.*ttest(allrecon1',0,'tail','right')';
% figure,ft_plot_mesh(source, 'vertexcolor', hbo_masked/max(abs(hbo_masked)), 'colormap', redblue);view(0,0);title('\DeltaHbO, p<0.05')
figure,ft_plot_mesh(source, 'vertexcolor', hbo_masked/max(abs(avg_recon1)), 'colormap', redblue);view(0,0);title('\DeltaHbO, p<0.05')
caxis([-0.8,0.8]); camlight headlight
saveas(gcf, 'Figures/HbO_masked.fig')
print Figures/HbO_masked.png -dpng -r300

figure,ft_plot_mesh(source, 'vertexcolor', avg_recon2/max(abs(avg_recon2)), 'colormap', redblue);view(0,0);title('\DeltaHb')
caxis([-0.8,0.8]); camlight headlight
saveas(gcf, 'Figures/Hb.fig')
print Figures/Hb.png -dpng -r300
hb_masked = avg_recon2.*ttest(allrecon2',0,'tail','left')';
% figure,ft_plot_mesh(source, 'vertexcolor', hb_masked/max(abs(hb_masked)), 'colormap', redblue);view(0,0);title('\DeltaHb, p<0.05')
figure,ft_plot_mesh(source, 'vertexcolor', hb_masked/max(abs(avg_recon2)), 'colormap', redblue);view(0,0);title('\DeltaHb, p<0.05')
caxis([-0.8,0.8]); camlight headlight
saveas(gcf, 'Figures/Hb_masked.fig')
print Figures/Hb_masked.png -dpng -r300
% 
% allrecon_hbt = allrecon1 - allrecon2;
% avg_recon_hbt = nanmean(allrecon_hbt, 2);
% figure,ft_plot_mesh(source, 'vertexcolor', avg_recon_hbt, 'colormap', redblue);view(0,0);title('\DeltaHbT')
% caxis([-max(abs(avg_recon_hbt)), max(abs(avg_recon_hbt))]); camlight headlight
% figure,ft_plot_mesh(source, 'vertexcolor', avg_recon_hbt.*ttest(allrecon_hbt')', 'colormap', redblue);view(0,0);title('\DeltaHbT, p<0.05')
% caxis([-max(abs(avg_recon_hbt)), max(abs(avg_recon_hbt))]); camlight headlight

% save('groupnirs_leftright.mat', 'allrecon1', 'allrecon2', 'avg_recon1', 'avg_recon2','hbo_masked', 'hb_masked', 'source', 'bad_subject')

% t_blk = tmp.t_blk;
% chanpos = tmp.chanpos;
% src_pos = tmp.src_pos;
% det_pos = tmp.det_pos;
% link = tmp.link;
% stim = tmp.stim;
% Fs = tmp.Fs;
% 
% stimvec = zeros(length(t_blk), 1);
% stimvec(25 + round((stim(1).data(1:20,1) - stim(1).data(1,1))*Fs)) = 1;
% stimvec(25 + round((stim(2).data(1:20,1) - stim(2).data(1,1))*Fs)) = 1;
% h = canonicalHRF(25, Fs);
% convolved = filter(h, 1, stimvec);
% designmat = [ones(length(t_blk),1), (1:length(t_blk))', convolved];
% % blk_avg(isnan(blk_avg)) = 0;
% tmp = pinv(designmat) * blk_avg;
% groupbeta = tmp(3,:)';
% A = jacobian_full;
% LL=spdiags(sqrt(sum(A.^2)'+0.1),0,size(A,2),size(A,2));
% LL_inv = inv(LL);
% At = A*LL_inv;
% alpha = 0.01;
% recon11 = LL_inv*(At'/alpha*groupbeta - At'*((inv(eye(nchannels*2) + At*At')*At)*(At'/alpha*groupbeta)));
% recon_groupbeta = griddata(nirs_mesh.nodes(:,1),nirs_mesh.nodes(:,2),nirs_mesh.nodes(:,3),recon11(1:length(nirs_mesh.nodes)),source.pos(:,1),source.pos(:,2),source.pos(:,3));
% figure,ft_plot_mesh(source, 'vertexcolor', recon_groupbeta(1:length(source.pos)), 'colormap', redblue);view(0,0)
% caxis([-max(abs(recon_groupbeta(1:length(recon_groupbeta)/2))), max(abs(recon_groupbeta(1:length(recon_groupbeta)/2)))]);
