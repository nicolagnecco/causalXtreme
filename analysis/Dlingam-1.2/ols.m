function [Bols, olsdisturbancestd, cols] = ols(X,k)
% OLS function written by Patrik O. Hoyer

dims = size(X,1);

%Permute the variables to the causal order
X = X(k,:);

% Remember to subract out the mean
Xm = mean(X,2);
X = X - Xm*ones(1,size(X,2));

% Calculate covariance matrix
C = (X*X')/size(X,2);

% Do QL decomposition on the inverse square root of C
[Q,L] = tridecomp(C^(-0.5),'ql');

% The estimated disturbance-stds are one over the abs of the diag of L
olsdisturbancestd = 1./diag(abs(L));

% Normalize rows of L to unit diagonal
L = L./(diag(L)*ones(1,dims));

% Calculate corresponding B
Bols = eye(dims)-L;

% Also calculate constants
cols = L*Xm;

% Permute back to original variable order
ik = iperm(k);
Bols = Bols(ik, ik);
olsdisturbancestd = olsdisturbancestd(ik);
cols = cols(ik);

return
