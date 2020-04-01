function Best_lasso = dolasso( X )

dims = size( X, 1 );
Best_lasso = zeros( dims, dims );

for i = 1 : dims
    y_lasso = X( i, : );
    theothers_lasso = 1:1:dims;
    theothers_lasso( i ) = [];
    X_lasso = X( theothers_lasso, : );
    
    % Run LASSO
    [b_lasso, fitinfo_lasso] = lasso(X_lasso',y_lasso');
    BICs = computeBIC( b_lasso, fitinfo_lasso, X_lasso, y_lasso );
    [~, bestBICindex] = min( BICs );
    b_lasso_out = b_lasso(:,bestBICindex)';
    
    Best_lasso( i, theothers_lasso ) = b_lasso_out;
end
