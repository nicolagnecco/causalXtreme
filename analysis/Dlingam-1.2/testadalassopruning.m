function testadalassopruning
% Version: 0.1
% Shohei Shimizu (9 Dec 2010)

% Set the randseed to the same each time
randseed = 0;
fprintf('Using randseed: %d\n',randseed);
rand('seed',randseed);
randn('seed',randseed);

p = 3
n = 100

nCV = 5

B = [ [ 0 0 0 ]; [ 1 0 0 ]; [ 0 1 0 ]]
A = inv( eye( p ) - B )

E = randn(p,n);
E = sign(E).*E.^(2);

X = A * E;

[Best, stdeest, ciest, kest] = Dlingam( X );
Best

Bpruned = adalassopruning(X, kest, nCV, Best )