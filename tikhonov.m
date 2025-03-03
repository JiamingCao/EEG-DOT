function [xhat, inv_op, res] = tikhonov(A, alpha, y)
% xhat = argmin_x ||Ax - y||^2 + ||x||^2
m = size(A,1);
n = size(A,2);
if m > n
    inv_op = inv(A'*A + alpha*eye(n))*A';
    if exist('y') && ~isempty(y)
        xhat = inv_op * y;
    else
        xhat = nan;
    end
    if nargout == 3
        res = norm(A*xhat - y);
    end
else
    M = A/sqrt(alpha);  
%   inv_op = (eye(n) - M'*(inv(eye(m) + M*M')*M))*(A'/alpha);
    inv_op = A' - M'*(inv(eye(m) + M*M')*(M*A'));
    inv_op = inv_op/alpha;
    if exist('y') && ~isempty(y)
        xhat = inv_op * y;
    else
        xhat = nan;
    end
    if nargout == 3
        res = norm(A*xhat - y);
    end
%     xhat = ((eye(n) - M'*(inv(eye(m) + M*M')*M))*(A'/alpha*y));
end
