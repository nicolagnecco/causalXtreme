function [B, stde, ci, K ] = pwDling(X, varargin)
% pwDling - performs pairwise~LiNGAM algorithm of Hyvarinen and Smith (2013,
% JMLR) in the DirectLiNGAM framework of Shimizu et al. (2011, JMLR).
%
% Before using this function, download Code for Pairwise Causality Measures from
% "http://www.cs.helsinki.fi/u/ahyvarin/code/pwcausal/,
% and add a path to "pwling" containing the extracted files.
%
% SYNTAX:
% [B, stde, ci, K] = pwDlin(X); No prior knowledge on the structure
% or
% [B, stde, ci, K] = pwDling(X, 'pk', Aknw); Some prior knowledge available
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
% Version: 0.9
% Shohei Shimizu (4 Oct 2015)
% Based on Dlingam code ver. 1.2

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
        index = findindex( X, candidates, U_K );
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

function index = findindex( X, candidates, U_K )

p = size(X,1);

% calculate T
T = NaN(1,p);

minT = -1; %% SS (24 Sep 2010)

for j = candidates
    
    if minT == -1 %% SS (24 Sep 2010) Input: minT, X, R, j
        T(j) = 0;
        for i = setdiff(U_K, j)
            LR = pwling( [ X(i,:); X(j,:) ], 1 );
            T(j) = T(j) + min( 0, LR(2,1) )^2;
        end
        minT = T(j);
    else
        T(j) = 0;
        for i = setdiff(U_K, j)
            LR = pwling( [ X(i,:); X(j,:) ], 1 );
            T(j) = T(j) + min( 0, LR(2,1) )^2;
            
            if T(j) > minT
                T(j) = Inf;
                break;
            end
            
        end
        minT = min( [ T(j), minT ] );
    end %% SS (24 Sep 2010) Output: minT, T(j)
    
end

% find argmin T
[~, index] = min(T);

end

