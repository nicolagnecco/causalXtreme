function b_ada_out = adalasso( y, X, nCV )
%
% Version: 0.99
% Shohei Shimizu (16 Jun 2012)

% Center variables
[ p, n ] = size( X );
y = y - mean( y, 2 ) * ones( 1, n );
X = X - mean( X, 2 ) * ones( 1, n );

%%%%%%%%%%%%% Do adaptive lasso %%%%%%%%%%%%%

% testgammas = [ 0.5, 1, 2 ];
testgammas = 1;
b_adaS = zeros( length( testgammas ), p );
MSES = zeros( length( testgammas ), 1 );

for i = 1 : length( testgammas )
    gamma = testgammas(i);
    
    % 1. Define X** (double-ast)
    b_ols = regress(y',X');
    w = 1 ./ ( abs( b_ols ).^(gamma) );
    
    Xdast = ( diag( w )^(-1) ) * X;
    
    % 2. Do lasso    
    [b fitinfo] = lasso(Xdast',y','CV',nCV,'Standardize','false');
    lam = fitinfo.IndexMinMSE; % find index of Minimum MSE
    b_opt = b(:,lam)';
    MSE = fitinfo.MSE( 1, lam );
    
    % 3. Compute b_ada
    
    b_ada = (( diag( w ) )^(-1) * b_opt')';
    
    % 4. Collect results
    
    b_adaS( i, : ) = b_ada;
    MSES( i, 1 ) = MSE;
    
end

[ ~, a ] = min( MSES );
b_ada_out = b_adaS( a,: );
