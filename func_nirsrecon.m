function func_nirsrecon(folder, source, jacobian_full, nirs_mesh)
file = dir([folder, '/*.snirf']);
rawdata = loadsnirf([folder, '/', file(end).name]);
Fs = 1/(rawdata.nirs.data.time(2) - rawdata.nirs.data.time(1)); % Hz

% Truncate data between first and last marker
t1 = max(1, rawdata.nirs.stim(1).data(1,1)-20);
t2 = min(rawdata.nirs.stim(1).data(end,1)+30, rawdata.nirs.data.time(end));
t1_idx = round(t1*Fs);
t2_idx = round(t2*Fs);
for i=1:length(rawdata.nirs.stim)
    rawdata.nirs.stim(i).data(:,1) = rawdata.nirs.stim(i).data(:,1) - t1;
end
time = rawdata.nirs.data.time(t1_idx:t2_idx) - rawdata.nirs.data.time(t1_idx);
rawdata.nirs.data.time = time;
rawdata.nirs.data.dataTimeSeries = rawdata.nirs.data.dataTimeSeries(t1_idx:t2_idx, :);

link = [extractfield(rawdata.nirs.data.measurementList,'sourceIndex');
    extractfield(rawdata.nirs.data.measurementList,'detectorIndex')]';
nchannels = size(link,1)/2;

time = rawdata.nirs.data.time;
stim = rawdata.nirs.stim;

% optode locations
det_pos = rawdata.nirs.probe.detectorPos2D;
src_pos = rawdata.nirs.probe.sourcePos2D;

% Convert to optical densities
nirs = rawdata.nirs.data.dataTimeSeries;
% nirs = TDDR(nirs, Fs);  % motion correction
% nirs = PCAFilter(nirs, 0.8);
dOD = bsxfun(@minus, log(nirs), log(median(nirs)));
dOD = TDDR(dOD, Fs);
% dOD = [TDDR(dOD(:,1:nchannels), Fs), TDDR(dOD(:,nchannels+1:end), Fs)];
% dOD = PCAFilter(dOD, 5);
% Detrend using 2nd order polynomial
for i=1:size(dOD, 2)
    polycoeff = polyfit([1:size(dOD, 1)]', dOD(:,i), 2);
    dOD(:,i) = dOD(:,i) - polyval(polycoeff, [1:size(dOD, 1)]');
end

win = 30;
% for j=1:size(dOD, 2)
%     tmp = dOD(:,j);
%     std_tmp = std(tmp);
%     for i=1:size(dOD,1)-win
%         if max(tmp(i:i+win))-min(tmp(i:i+win))>3*std_tmp && abs(tmp(i)-tmp(i+win))<std_tmp
%             tmp(i:i+win)=tmp(i)+[0:win]*(tmp(i+win)-tmp(i))/win;
%         end
%     end
%     dOD(:,j) = tmp;
% end

% for j=1:size(dOD, 2)
%     tmp = dOD(:,j);
%     for i=1:size(dOD,1)-win
%         if max(tmp(i:i+win))-min(tmp(i:i+win))>0.1 && abs(tmp(i)-tmp(i+win))<0.05
%             tmp(i:i+win)=tmp(i)+[0:win]*(tmp(i+win)-tmp(i))/win;
%         end
%     end
%     dOD(:,j) = tmp;
% end

warning('off', 'signal:findpeaks:largeMinPeakHeight');
for j=1:size(dOD,2)
    tmp = dOD(:,j);
    sqdiff = zeros(size(dOD,1)-win,1);
    for i=1:size(dOD,1)-win
        sqdiff(i) = sum(diff(tmp(i:i+win)).^2);
    end
    [~,loc]=findpeaks(sqdiff,'MinPeakDistance',30,'MinPeakHeight',0.02);
    for i=1:length(loc)
        tmp(loc(i):loc(i)+win)=tmp(loc(i))+[0:win]*(tmp(loc(i)+win)-tmp(loc(i)))/win;
    end
    dOD(:,j) = tmp;
end

% filter
[b,a] = butter(3, 0.01/(Fs/2), 'high');
dOD_highpass = filtfilt(b, a, dOD);
% [b,a] = butter(3, [0.01, 0.1]/(Fs/2));
% dOD_filt = filtfilt(b, a, dOD);

% bad = [];
% for i=1:nchannels
%     if std(dOD_highpass(:,i)) > 0.08 || std(dOD_highpass(:,i+nchannels)) > 0.08
%         bad = [bad; i];
%     end
% end
[freq,amp] = get_spec(Fs, mean(dOD_highpass(:,std(dOD_highpass)<0.08), 2));
amp(freq<0.8 | freq>1.5) = 0;
[~, heartrate_idx] = max(amp);
heartrate = freq(heartrate_idx);
fprintf('Heart rate: %fHz\n', heartrate);
snr = zeros(nchannels, 1);
for i=1:nchannels
    [freq,amp] = get_spec(Fs, dOD_highpass(:,nchannels+i));
    signal_pow = sum(amp(freq>heartrate-0.1 & freq<heartrate+0.1));
    noise_pow = sum(amp(freq>1.5*heartrate-0.1 & freq<1.5*heartrate+0.1));
    snr(i) = signal_pow / noise_pow;    % estimate SNR based on heart rate peak
end
bad = find(snr<1.5);
fprintf('Bad channels: %d\n', length(bad));
if length(bad)>nchannels/2
    warning('High number of noise channels')
end

%% GLM
stimvec = zeros(length(time), 1);
stimvec(round(stim(1).data(:,1)*Fs)) = 1;
% stimvec(round(stim(2).data(:,1)*Fs)) = 1;
h = canonicalHRF(32, Fs);
convolved = filter(h, 1, stimvec);

if ~isfield(rawdata.nirs, 'aux')
    designmat = [(1:length(time))', ones(length(time),1), convolved, mean(dOD_highpass, 2)];
%     designmat = [legendre(time, 3), convolved, mean(dOD_highpass, 2)];
    fprintf('No accelorometer found\n')
else
    [b,a] = butter(3, 0.01/(Fs/2), 'high');
    tmp = split([folder, '/', file(end).name], '.');
    load([tmp{1}, '.nirs'], '-mat', 'aux'); % Accelorometer data is at higher sampling frequency in .snirf, but correct in .nirs 
    accelorometer = [squeeze(aux(:,1,:)), squeeze(aux(:,2,:))];
    designmat = [(1:length(time))', ones(length(time),1), convolved, filtfilt(b,a,accelorometer(t1_idx:t2_idx, :))];
%     designmat = [legendre(time, 3), convolved, filtfilt(b,a,accelorometer(t1_idx:t2_idx, :))];
end
% tmp = pinv(designmat) * dOD_filt;

Pmax = ceil(4*Fs);
[tmp,stats] = glm_ar_irls(dOD_highpass-mean(dOD_highpass), designmat, Pmax);    % GLM after prewhitening; functions from Ted Huppert toolbox
% [tmp,stats] = glm_ar_irls(dOD-mean(dOD), designmat, Pmax);
% good = setdiff(1:nchannels, bad);
% tmp = glm_ar_irls(dOD_highpass(:,[good,good+nchannels])-mean(dOD_highpass(:,[good,good+nchannels])), designmat, Pmax);
beta = tmp(3,:)';
% beta = tmp(5,:)';

%% Reconstruction
A = jacobian_full;
A([bad;bad+nchannels], :) = [];
beta([bad; bad+nchannels]) = [];
LL=spdiags(sqrt(sum(A.^2)'+0.1),0,size(A,2),size(A,2)); % spatially varying regularization
LL_inv = inv(LL);
At = A*LL_inv;
alpha = 0.01;
% Acceleration using Woodbury equality
recon11 = LL_inv*(At'/alpha*beta - At'*((inv(eye(size(A,1)) + At*At')*At)*(At'/alpha*beta)));
F = scatteredInterpolant(nirs_mesh.nodes(:,1),nirs_mesh.nodes(:,2),nirs_mesh.nodes(:,3),recon11(1:length(nirs_mesh.nodes)));
recon1 = F(source.pos(:,1),source.pos(:,2),source.pos(:,3));
F = scatteredInterpolant(nirs_mesh.nodes(:,1),nirs_mesh.nodes(:,2),nirs_mesh.nodes(:,3),recon11(length(nirs_mesh.nodes)+1:end));
recon2 = F(source.pos(:,1),source.pos(:,2),source.pos(:,3));
figure,ft_plot_mesh(source, 'vertexcolor', recon1, 'colormap', redblue);view(0,0);camlight headlight
caxis([-max(abs(recon1)), max(abs(recon1))]);
saveas(gcf, [folder, '/recon4.fig']);
saveas(gcf, [folder, '/recon4.png']);
close gcf
figure,ft_plot_mesh(source, 'vertexcolor', recon2, 'colormap', redblue);view(0,0);camlight headlight
caxis([-max(abs(recon2)), max(abs(recon2))]);
saveas(gcf, [folder, '/recon_hb4.fig']);
saveas(gcf, [folder, '/recon_hb4.png']);
close gcf

save([folder, '/processed_recon4'],'bad','beta','designmat','heartrate','stats','det_pos','Fs','link','nchannels','recon1','recon2','src_pos','stim','time')
fprintf('Done processing\n')

