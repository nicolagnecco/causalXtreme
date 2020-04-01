function Aest = estA(X,k)
% Estimate A without computing matrix inversion
% Also A can be computed by inv(I-B).
% 
% Version: 0.99
% Shohei Shimizu (9 Dec 2010)

p = size(X,1);
Aest = eye( p );

%Permute the variables to the causal order
X = X(k,:);

% Remember to subract out the mean
Xm = mean(X,2);
X = X - Xm*ones(1,size(X,2));

% Do regression analysis
R = X;

for j = 1 : p
    
%     Aest
    
    Cov = cov(R');
    
    for i = j : p
        if i ~= j
            Aest(i,j) = Cov(i,j)/Cov(j,j);
            R(i,:) = R(i,:) - Cov(i,j)./Cov(j,j)*R(j,:);
        end
    end
    
end

% Permute back to original variable order
ik = iperm(k);
Aest = Aest(ik, ik);

return;