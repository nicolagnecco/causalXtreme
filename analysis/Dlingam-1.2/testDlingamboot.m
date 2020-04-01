function testDlingamboot
% Version: 0.1
% Shohei Shimizu (9 Dec 2010)

% Set the randseed to the same each time
randseed = 0;
fprintf('Using randseed: %d\n',randseed);
rand('seed',randseed);
randn('seed',randseed);

% parameters
nboot = 1000
threshold = 0.05

% Generate data
p = 3
n = 100

B = [ [ 0 0 0 ]; [ 1 0 0 ]; [ 0 1 0 ]]
A = inv( eye( p ) - B )

E = randn( p, n );
E = sign( E ) .* E.^2;

X = A * E;

[Best, stdeest, ciest, kest] = Dlingam(X)
Aest = estA(X,kest) % or inv(eye(p)-Best);

[ Bsig, Asig ] = Dlingamboot( X, nboot, threshold )
% [ Bsig, Asig, LB, UB, LA, UA ] = Dlingamboot( X, nboot, threshold ) % Also outputs the lower bounds and upper bounds