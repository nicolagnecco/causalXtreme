function BICs = computeBIC( b, fitinfo, X, y )
% Written by Shohei Shimizu and Kittitat Thamvitayakul

L = length( fitinfo.Lambda ); % Number of Lambda parameters

% Center variables
[ p, n ] = size( X );
y = y - mean( y, 2 ) * ones( 1, n );
X = X - mean( X, 2 ) * ones( 1, n );

% Compute BICs
if fitinfo.Alpha == 1 % lasso
    
    sigmae2 = sum( ( y' - X' * pinv( X' ) * y' ).^2 ) / n;
    vdf = sum( b' ~= 0 , 2 );
    vy = repmat(y,L,1);
    BICs = ( sum ( ( vy - b' * X ).^2 , 2 ) + log( n ) * sigmae2 * vdf )';
    
else % Elastic net
    
    vdf = zeros( L, 1 );
    deltaS = fitinfo.Lambda * ( 1 - fitinfo.Alpha ) / 2;
    
    beta_ridge_tmpS = ridge( y', X', deltaS, 0 ); % p+1 x L
    beta_ridgeS = beta_ridge_tmpS( 2 : p+1, : ); % p x L
    vy = repmat(y,L,1); % L x n
    sigmae2 = sum( ( vy - beta_ridgeS' * X ).^2, 2 ) / n; % L x 1
    
    for i = 1 : L
        
        delta = deltaS( i );
        
        X_A = X(  b( :, i ) ~= 0 , : );
        I_A = eye( size( X_A, 1 ) );
        
        vdf( i ) = trace( X_A' / ( X_A * X_A' + delta * I_A ) * X_A );
                
    end
    
    BICs = sum( ( vy - b' * X ).^2, 2 ) + log( n ) * sigmae2 .* vdf; % L x 1
    BICs = BICs';
 
end

