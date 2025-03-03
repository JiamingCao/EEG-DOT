function data = PCAFilter( data, ncomp)

% resample data
d = data;

% remove mean
m = mean(d,1);
d = bsxfun(@minus, d, m);

% svd
[u, s, v] = svd(d,'econ');
s = diag(s);

% remove n components
s(1:ncomp) = 0;
d = u*diag(s)*v';

% add mean back
d = bsxfun(@plus, d, m);

% put back
data = d;
end

