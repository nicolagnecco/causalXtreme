function [B, stde, ci, K] = Dlingam(X, varargin)
% Dlingam - perform the DirectLiNGAM algorithm for the paper "A direct
% method for learning a linear non-Gaussian structural equation model" by
% Shimizu et al.
%
% Before using this function, download KernelICA codes from
% "http://www.di.ens.fr/~fbach/", extract the files,
% and add a path to "kernel-ica1_2" containing the extracted files.
%
% SYNTAX:
% [B, stde, ci, K] = Dlingam(X); No prior knowledge on the structure
% or
% [B, stde, ci, K] = Dlingam(X, 'pk', Aknw); Some prior knowledge available
%
% INPUT:
% X     - Data matrix: each row is an observed variable, each
%         column one observed sample. The number of columns should
%         be far larger than the number of rows.
% Aknw     - Matrix of prior knowledge
%  Aknw(i,j) is 0 if x_j does NOT have a directed path to x_i
%  Aknw(i,j) is 1 if x_j has a directed path to x_i
%  Aknw(i,j) is -1 if no prior knowledge available
% (A directed path from x_j to x_i is a sequence of directed edges such that
% x_i is reachable from x_j.)
%
% OUTPUT:
% B     - Matrix of estimated connection strenghts
% stde  - Standard deviations of disturbance variables
% ci    - Constants
% K     - An estimated causal ordering
%
% Version: 0.1
% Takanori Inazumi (21 Jun 2010)
% Version: 1.0
% Updated by Shohei Shimizu (24 Sep 2010, 8 Nov 2010, 10 Dec 2010)
% We thank Yasuhiro Sogawa

% 1. %
Xorig = X;

[p, n] = size(X); % p: N of var, n: N of data points

% Center variables
X = X - mean(X,2)*ones(1,n);

% Prior knowledge matrix Aknw % SS (24 Sep 2010)
M = -1 * ones( p, p );
for i = 1:2:length(varargin) - 1
    switch varargin{i}
        case 'pk' % Prior Knowledge
            Aknw = varargin{i + 1};
            if sum(sum( Aknw~=0 & Aknw~=1 & Aknw~=-1)) ~= 0
                error('Aknw includes improper elements.');
            end
            M = Aknw;
        otherwise
            error('unknown argument name [%s]', varargin{i});
    end
end

% erase diagonal elements of M
M(logical(eye(p,p))) = 0;

K = zeros(1,p);
U_K = 1:p; % U-K

% 2. %D
for m = 1:p-1
    % find exogenous by using M
    exogenous = find(sum(M' == 0) == p-m+1); % length(U_K) == p-m+1
    if isempty(exogenous)
        endgenous = find(sum(M' == 1) > 0);
        candidates = setdiff(U_K, endgenous);
    else
        candidates = exogenous;
    end
    
    % (a) %
    % calculate R
    R = computeR( X, candidates, U_K, M );
    
    % skip exogenous finding if it is found
    if length(candidates) == 1
        index = candidates;
    else
        % find exogenous
        index = findindex( X, R, candidates, U_K );
    end
    
    % (b) %
    K(m) = index;
    U_K(U_K == index) = [];
    M(:,index) = NaN;
    M(index,:) = NaN;
    
    % (c) %
    X = R(:,:,index);
end

% 3. %
K(p) = U_K;

% 4. %
[B, stde, ci] = ols(Xorig, K);

end

% function y = g(x)
%
% y = tanh(x);
%
% end

function [J] = my_call_contrast(x)
% Author: Yasuhiro Sogawa
% Modified by SS (27 Sep 2010)
% my_call_contrast - set parameters of KernelICA.
% and call a contrast function employed in KernelICA.
% The details of the parameters are shown in section 4.5,
% "Kernel Independent Component Analysis" (F. R. Bachand and M.I.Jordan).

[m,N]=size(x);

% set the parameters
contrast='kgv';
% contrast='kcca';
if N < 1000
    sigma=1;
    kappa=2e-2;
else % Added by SS (24 Sep 2010)
    sigma = 1/2;
    kappa = 2e-3;
end

kernel='gaussian';

mc=m;
kparam.kappas=kappa*ones(1,mc);
kparam.etas=kappa*1e-2*ones(1,mc);
kparam.neigs=N*ones(1,mc);
kparam.nchols=N*ones(1,mc);
kparam.kernel=kernel;
kparam.sigmas=sigma*ones(1,mc);

% Commented out by SS (24 Sep 2010)
% % scales data
% covmatrix=x*x'/N;
% sqrcovmatrix=sqrtm(covmatrix);
% invsqrcovmatrix=inv(sqrcovmatrix);
% x=invsqrcovmatrix*x;

% perform contrast function
J = contrast_ica(contrast,x,kparam);

end

function R = computeR( X, candidates, U_K, M )

[p,n] = size( X );
R = zeros(p,n,p);
Cov = cov(X');

for j = candidates
    for i = setdiff(U_K, j)
        % skip residue calculation by using M
        if M(i,j) == 0
            R(i,:,j) = X(i,:);
        else
            R(i,:,j) = X(i,:) - Cov(i,j)/Cov(j,j)*X(j,:);
        end
    end
end

end

function index = findindex( X, R, candidates, U_K )

p = size(X,1);

% calculate T
T = NaN(1,p);

minT = -1; %% SS (24 Sep 2010)

for j = candidates
    % Commentted out by SS (24 Sep 2010)
    %             T(j) = 0;
    %             for i = setdiff(U_K, j)
    %                 Corr = corrcoef([g(R(i,:,j)); X(j,:); R(i,:,j); g(X(j,:))]');
    %                 T(j) = T(j) + abs(Corr(1,2)) + abs(Corr(3,4));
    %             end
    
    if minT == -1 %% SS (24 Sep 2010) Input: minT, X, R, j
        T(j) = 0;
        for i = setdiff(U_K, j)
            J = my_call_contrast([R(i,:,j); X(j,:)]); %using kernel based independence measure
            T(j) = T(j) + J;
        end
        minT = T(j);
    else
        T(j) = 0;
        for i = setdiff(U_K, j)
            J = my_call_contrast([R(i,:,j); X(j,:)]); %using kernel based independence measure
            T(j) = T(j) + J;
            
            if T(j) > minT
                T(j) = Inf;
                break;
            end
            
        end
        minT = min( [ T(j), minT ] );
    end %% SS (24 Sep 2010) Output: minT, T(j)
    
end

% find argmin T
[minval, index] = min(T);

end
