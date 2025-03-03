function [lambda,Beta,inv_op]=REML(Y,X,Beta_prior,Qn,Qp,maxIter,lambda,jump, print_flag)
% Restricted Maximum Likelihood code.  This program was modified from the
% SPM function spm_reml to be optimized for the dimensions of optical data.
% For reference see
%       Friston KJ. Statistical parametric mapping : the analysis of functional
%       brain images. London: Academic; 2007.
%
% Copied from Ted Huppert's toolbox, with modification

if(nargin<3 || isempty(Beta_prior))
    Beta_prior=zeros(size(X,2),1);
end
if(nargin<4 || isempty(Qn))
    Qn{1}=speye(size(Y,1),size(Y,1));
end
if(nargin<5 || isempty(Qp))
    Qp{1}=speye(size(X,2),size(X,2));
end
if(nargin<6 || isempty(maxIter))
    maxIter=50;
end
if ~exist('print_flag') || isempty(print_flag)
    print_flag = true;
end

% [~,s,~]=nirs.math.mysvd(Y);
% s1=max(diag(s));
% Y=Y/s1;
% [~,s,~]=nirs.math.mysvd(X);
% s2=max(diag(s));
% X=X/s2;
% 
% Beta_prior=Beta_prior*s2/s1;

Beta=[];
Stats=[];

if(exist('jump')~=1)
    jump=false;
end
    
Beta_priorO=Beta_prior;

[m,n] = size(X);
[U,S,V]= mysvd(full(X));
s = diag(S);
tol = max(m,n) * eps(max(s));
r = sum(s > tol);
s = diag(s(1:r));
X = U(:,1:r)*s*V(:,1:r)'; 

tolr=eps(full(max(X(:))))*max(size(X))*10;


if(~exist('maxIter'))
    maxIter=150;  %Max # of iterations of REML code
end


if(size(X,1)<size(X,2))
   % If X is not full rank, then we can make this MUCH faster
   %  noting that the hyper-parameters of the reduced problem are the same as
   %  the orginal forward model.  Note- this was not done in the paper, but is
   %  a worthwhile extension in future work.  This trick was discovered after
   %  we finished the paper in prep of this demo.

    %  Y = U*S*V'*Beta
    %  Y = U*(S*V'*Beta) --> Y=U*S*Beta2;
    %  cov(Beta2) = S*V'*Q*V*S';
    %[V,S,U]=svd(full(X'),0);
    [U,S,V]=mysvd(full(X));
    
    for idx=1:length(Qp)
        Qp2{idx}=V'*Qp{idx}*V;
    end
    Beta_prior=V'*Beta_prior;
    [lambda,Beta,~]=REML(Y,U*S,Beta_prior,Qn,Qp2,maxIter);

    if(jump | nargout==1)
        return
    end
    
    %lambda is right, but the Stats are not directly related to the ones we want.  So we
    %recompute. 
%     Cn=tolr*speye(size(Qn{1},1));
%     for idx=1:length(Qn)
%         Cn=Cn+Qn{idx}*exp(lambda(idx));
%     end
%     Cp=tolr*speye(size(Qp{1},1));
%     for idx2=1:length(Qp)
%         Cp=Cp+Qp{idx2}*exp(lambda(idx+idx2));
%     end
%     Ce=blkdiag(Cn,Cp);
%     Ce=Ce+speye(size(Ce))*tolr;

%     iCn=inv((Cn+speye(size(Cn))*tolr));
%     iCp=inv((Cp+speye(size(Cp))*tolr));
%     iCe = blkdiag(iCn,iCp);
%     X2 = [X; speye(size(X,2))];
    
%     R=speye(size(X2,1))-X2*pinv(full(X2'*iCe*X2+10*eps(1)*speye(size(X2,2),size(X2,2))))*X2'*iCe;
        
else

    %%Else, run the normal model
    Y2=Y;
    X2=X;
    
    %% Set up heirarchical model
    Y = [Y; sparse(size(X,2),size(Y,2))];
    X = [X; speye(size(X,2))];

    
    %Set up the extended covariance model by concatinating the measurement
    %and parameter noise terms
    Q=cell(length(Qn)+length(Qp),1);
    for idx=1:length(Qn)
        Q{idx}=blkdiag(Qn{idx},sparse(size(Qp{1},1),size(Qp{1},2))); % Build block diagonal matrix from Qn & Qp matrices
    end
    for idx2=1:length(Qp)
        Q{idx+idx2}=blkdiag(sparse(size(Qn{1},1),size(Qn{1},2)),Qp{idx2});
    end

%     for idx=1:length(Q)
%         Q{idx}=Q{idx}.*(abs(Q{idx})>max(abs(Q{idx}(:)))/1E4);
%     end
    
    
    if(~exist('lambda') || isempty(lambda))
        lambda = ones(length(Q),1);  %Initial guess of lambda. 
    end
    tol=1E-4;  %REML goes till tolorence (or max iter) 
    
    t=32;
    dF    = Inf;  %Initial decent 
    cnt=0;  %This is a bookkeeping param for display purposes

    %% This function was adapted from the SPM function spm_reml.m
    %  which is part of the SPM toolkit (http://www.fil.ion.ucl.ac.uk/spm)
    %

    for iter=1:maxIter
        Ce = tolr*speye(size(Q{1},1));  %Make sure it stays in numerical precision

        % E-step:

        %Bound lambda to avoid numerircal prec. issues
        lambda=max(lambda,log(tolr));
        lambda=min(lambda,log(1/tolr));

        for i = 1:length(Q)
            Ce = Ce + Q{i}*exp(lambda(i));
        end

  %      iCe=blkdiag(inv(Ce(1:end/2,1:end/2)),inv(Ce(1+end/2:end,1+end/2:end)));
        iCe = pinv(full(Ce));
        Xt_iCe = X' * iCe;
        Xt_iCe_X = Xt_iCe * X;
        C_beta_y = pinv(full(Xt_iCe_X));  %Estimate of covariance of beta given the measurements

        
        % M-step:
        P = iCe - (iCe*X)*C_beta_y*Xt_iCe;
        PY=P*Y;
        for i=1:size(Q,1)
            PQ_i{i}=P*Q{i};
        end

        for i = 1:size(Q,1)
            PQ = PQ_i{i};
            PQt=PQ';
            [u,s]=mysvd(Q{i});
            nn= 0.5*norm(PY'*u*sqrt(s))^2;
            g(i,1) = -0.5*trace(PQ)*exp(lambda(i)) +nn*exp(lambda(i));
            for j = i:size(Q,1)
                PQj = PQ_i{j};
                H(i,j) = -0.5*sum(sum(PQt.*PQj))*exp(lambda(i)+lambda(j));
                H(j,i)=H(i,j);
            end
        end

        %Now update the lambda.  dLambda = -inv(H)*g
        I=eye(size(H,1));
%        warning('off','MATLAB:nearlySingularMatrix');
        
       dL = (expm(H*t) - I)*pinv(full(H)+eye(size(H))*1E-10)*g;
       %dL= pinv(H)*g;
       lambda = lambda + dL;
        
        if(any(isnan(dL)))
            lambda=ones(size(lambda));
            disp('starting over');
            continue;
        end

        df    = g'*dL;
        if df > dF - exp(-2.5), t = max(2,t/2); end %retune the regularization if req., originally exp(-4)
        dF    = df;
        
        if print_flag
        for c=1:cnt, fprintf('\b'); end
        cnt=fprintf('%-5s: %i %5s%e','  ReML Iteration',iter,'...',full(dF));
        end

        if dF < tol && iter>5, break; end

    end
    
%     R=speye(size(X,1))-X*pinv(full(X'*iCe*X))*X'*iCe;

    
    X=X2;
    Y=Y2;
    
end

lambda=max(lambda,log(tolr));
lambda=min(lambda,log(1/tolr));

if(jump)
        return
end
clear U V X2 iCe iCp Qp2 iCe R S Cn Cp iCn Beta

%[Beta,stdx,mse]=regress_wLS(Y,X,Beta_priorO,Qn,Qp,lambda);
%[Beta,stdx,mse]=invl(Y,X,Qn,Qp,lambda,true);    


Cn = tolr*speye(size(Qn{1},1));  %Make sure it stays in numerical precision
for i = 1:length(Qn)
    Cn = Cn + Qn{i}*exp(lambda(i));
end
Cp = tolr*speye(size(Qp{1},1));  %Make sure it stays in numerical precision
for i = 1:length(Qp)
    Cp = Cp + Qp{i}*exp(lambda(i+length(Qn)));
end
% iCn=pinv(full(Cn));
% iCp=pinv(full(Cp));
% XtXi = pinv(X'*iCn*X+iCp);
% Beta = XtXi*X'*iCn*Y;

% Speed up using Woodbury's
iCn=pinv(full(Cn));
tmp = pinv(Cn+X*Cp*X');
inv_op = Cp*X'*iCn - (Cp*X'*tmp)*(X*Cp*X'*iCn);
Beta = inv_op*Y;

%Now, put the final pieces together
%Beta= C_beta_y * Xt_iCe * Y;


% clear Beta_prior CY C_beta_y Cn Cp Qn Qp X X2 iCe residuals ybar yhat 

fprintf('\n');

return