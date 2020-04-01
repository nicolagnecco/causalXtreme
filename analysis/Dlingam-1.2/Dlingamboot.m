function [Bsig, Asig, LB, UB, LA, UA] = Dlingamboot(X, nboot, threshold)
% Computes bootstrap confidence intervals using the percentile method and return significant effects.
% [Bsig, Asig] = Dlingamboot(X, nboot)
% X: data matrix
% nboot: number of bootstrap replicates, say 1000 or 10000.
% threshold: significance level, say 0.05
% Bsig: Matrix that collects significant directed
% edges (direct causal effects) bij of B
% Asig: Matrix that collects significant total causal effects of aij of A =
% I-B
% 
% Version: 0.99
% Shohei Shimizu (8 Dec 2010)

% Set the randseed to the same each time
randseed = 0;
fprintf('Using randseed: %d\n',randseed);
rand('seed',randseed);
randn('seed',randseed);

% Do DirectLiNGAM
Best = Dlingam(X);
Aest = inv( eye( size( X, 1 ) ) - Best );

%%% Do bootstrapping %%%
% p: number of variables. n: sample size
[p,n] = size(X);

% Get bootstrap estimates
fprintf('Bootstrapping...\n');

% This collects bootstrap estimates
Bboot = zeros(p,p,nboot);
Aboot = zeros(p,p,nboot);

% Do bootstrapping
for booti = 1 : nboot
    
    if rem(booti,ceil(10)) == 0
        fprintf('[%1.0f]',booti);
    end
    
    % Generate bootstrap samples
    bootindex = ceil( rand(1,n) * n );
    Xboot = X(:,bootindex);
    
    % Compute B
    [ Bestboot, dummy1, dummy2, kestboot ] = Dlingam(Xboot);
    
    % Collect estimates
    Bboot(:,:,booti) = Bestboot;
    Aestboot = estA(Xboot, kestboot);% or Aestboot = inv( eye( p ) - Bestboot );
    Aboot(:,:,booti) = Aestboot;
    
end
fprintf('\n');

% Compute bootstrap confidence intervals
Bboots = sort(Bboot,3);
LB = Bboots(:,:,ceil(nboot*threshold/2));% Lower bound of B
UB = Bboots(:,:,ceil(nboot*(1-threshold/2)));% Upper bound of B

Aboots = sort(Aboot,3);
LA = Aboots(:,:,ceil(nboot*threshold/2));% Lower bound of A
UA = Aboots(:,:,ceil(nboot*(1-threshold/2)));% Upper bound of A

% Set bij whose confidence intervals contain zero to zero
Bsig = Best;
Bsig( LB <=0 & UB >=0 )=0;

% Set aij whose confidence intervals contain zero to zero
Asig = Aest;
Asig( LA <=0 & UA >=0 )=0;

return;


    