function h = canonicalHRF(duration, Fs, peakTime, uShootTime, peakDisp, uShootDisp, ratio)
% Canonical double-gamma hrf. Adapted from nirs-toolbox
if nargin<3
    peakTime    = 6;    % SPM12 standard parameters
    uShootTime  = 16;
    peakDisp    = 1;
    uShootDisp  = 1;
    ratio       = 1/6;
end
t = 0:1/Fs:duration;
a1 = peakTime;
a2 = uShootTime;
b1 = peakDisp;
b2 = uShootDisp;
c  = ratio;
h = b1^a1*t.^(a1-1).*exp(-b1*t)/gamma(a1) - c*b2^a2*t.^(a2-1).*exp(-b2*t)/gamma(a2);
h = h / sum(h);
