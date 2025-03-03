function x = TDDR(x,fs,maxN,T)
%%
%   Fishburn F.A., Ludlum R.S., Vaidya C.J., & Medvedev A.V. (2019).
%   Temporal Derivative Distribution Repair (TDDR): A motion correction
%   method for fNIRS. NeuroImage, 184, 171-179.
%   https://doi.org/10.1016/j.neuroimage.2018.09.025

% x inoput signal x = [Time, Channel]
% fs sampling frequency
% maxN is the maximum number of iterations to fit the noise
% T is the threshold for the interative estimation of the noise.
% 
% Adapted from Ted Huppert's nirs-toolbox
%

if nargin < 4
    T=[];
end
if nargin < 3
    maxN = [];
end

if isempty(maxN)
    maxN=100;
end
if isempty(T)
    T=eps(1);
end


DC = median(x);
x=x-DC;

%% Preprocess: Separate high and low frequencies
filter_cutoff = .5;
filter_order = 3;
Fc = filter_cutoff * 2/fs;
if Fc<1
    [fb,fa] = butter(filter_order,Fc);
    x_low = filtfilt(fb,fa,x);
else
    x_low = x;
end
x_high = x - x_low;

y = diff(x_low);
w = ones(size(y));
mu = inf;
for ii=1:maxN
    mu_post = mu;
    
    % compute current iteration
    mu = sum(w.*y)./sum(w);
    r = abs(y-mu);
    s = 1.4826.*median(r);
    d = r./(4.685*s);
    % if median(r) is for whatever reason 0:
%     d(isinf(d)) = 1;
%     d(isnan(d)) = 1;
    %prepare next iteration
    w = ((1 - d.^2).*(d < 1)).^2;
    
    if abs(mu_post-mu) < T
        disp('threshold')
        break
    end
end
if ii==maxN; disp('iteration');end
yd = w.*(y-mu);
x_corrected = cumsum([zeros([1,size(yd,2)]);yd]);
x_corrected = x_corrected - mean(x_corrected);
x = x_corrected + x_high;
x=x+DC;

end
