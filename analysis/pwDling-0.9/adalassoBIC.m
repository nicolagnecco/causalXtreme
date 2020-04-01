function b_ada = adalassoBIC( y, X, varargin )
%
% Version: 0.99
% Shohei Shimizu (16 Jun 2012)

ridgeornot = 0;
for i = 1:2:length(varargin) - 1
    switch varargin{i}
        case 'ridge'
            ridgeornot = varargin{i + 1};
        otherwise
            error('unknown argument name [%s]', varargin{i});
    end
end

% Center variables
n = size( X, 2 );
y = y - mean( y, 2 ) * ones( 1, n );
X = X - mean( X, 2 ) * ones( 1, n );

%%%%%%%%%%%%% Do adaptive lasso %%%%%%%%%%%%%

gamma = 1; % 0.5, 1, or 2

% 1. Define X** (double-ast)
if ridgeornot == 1
    b_ridge = doridge( y, X );
    b_ridge = b_ridge';
    w = 1 ./ ( abs( b_ridge ).^(gamma) );
else
    b_ols = regress(y',X');
    w = 1 ./ ( abs( b_ols ).^(gamma) );
end

Xdast = ( diag( w )^(-1) ) * X;

% 2. Do lasso
[b, fitinfo] = lasso(Xdast',y','Standardize',false);
BICs = computeBIC( b, fitinfo, Xdast, y );
[~, bestBICindex] = min( BICs );
b_opt = b(:,bestBICindex)';

% 3. Compute b_ada
b_ada = (( diag( w ) )^(-1) * b_opt')';

