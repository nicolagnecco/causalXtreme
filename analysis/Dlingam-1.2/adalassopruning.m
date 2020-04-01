function Bpruned = adalassopruning(X, k, nCV, B )
% Bpruned = adalassopruning(X, k, nCV, B )
% X: data matrix
% k: an estimated ordering
% nCV: nCV-fold cross-validation, say 5
% B: an estimated connection strength matrix
% 
% Version: 0.9
% Shohei Shimizu (8 Dec 2010)


% Get the dimension of X
dim = size( X, 1 );

% Permute the variables to the causal order k
Xp = X(k,:);
Bp = B(k,k);

% Do adaptive lasso

Badap = zeros( size( Bp ) );

for i = 1 : dim        
    if sum( Bp( i, : ) ~= 0 ) > 0

        yp = Xp( i, : ); 
        
        predictorsp = find( Bp( i, : ) ~= 0 );        
        
        Xpredp = Xp( predictorsp, : );        

        badap = adalasso( yp, Xpredp, nCV );
        
        Badap( i, 1:length(badap) ) = badap;

    end
end


% Permute back to original variable order
ik = iperm(k);
Bada = Badap(ik, ik);

Bpruned = B;
Bpruned( Bada == 0 ) = 0;

