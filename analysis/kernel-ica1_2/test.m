clc
rng(42)
n = 1e5;
X1 = randn(1, n);
X2 = X1 + randn(1, n);
Cov = cov(X1, X2);
R2 = X2 - Cov(1, 2) / Cov(1, 1) * X1;

x = [R2; X1];


x = x(:, 1:100);

% Call Kernel ICA
[m,N]=size(x);


% set the parameters
contrast='kgv';
% contrast='kcca';
if N < 1000
    sigma=1;
    kappa=2e-2;
else % Added by SS (24 Sep 2010)
    sigma = 1/2;
    kappa = 2e-3;
end

kernel='gaussian';

mc=m;
kparam.kappas=kappa*ones(1,mc);
kparam.etas=kappa*1e-2*ones(1,mc);
kparam.neigs=N*ones(1,mc);
kparam.nchols=N*ones(1,mc);
kparam.kernel=kernel;
kparam.sigmas=sigma*ones(1,mc);

tic
J = contrast_ica(contrast, x, kparam);
toc

%%
% csvwrite('../X_mat.csv', x')
%% function contrast_ica
clc
format long g
N=size(x,2); 		% number of data points
m=size(x,1);      % number of components
kappas=kparam.kappas;
etas=kparam.etas;
Rkappa=[];
sizes=[];
i = 1;
tic
[G,Pvec] =chol_gauss(x(i,1:100)/kparam.sigmas(i),1,N*kparam.etas(i));
toc

%% 
x = [1:30];
tic
[G, pvec] = chol_gauss(x, 1, 1);
toc

%%
G = [1:3; 3:5];
%% regularization (see paper for details)
[A,D]=eig(G'*G);
D=diag(D);

%%
indexes=find(D>=N*etas(i) & isreal(D)); %removes small eigenvalues
[newinds,order]=sort(D(indexes));

%%
order=flipud(order);

%%
neig=length(indexes);
indexes=indexes(order(1:neig));
if (isempty(indexes)), indexes=[1]; end
D=D(indexes);
V=G*(A(:,indexes)*diag(sqrt(1./(D))));
Us{i}=V;
Lambdas{i}=D;
Dr=D;
for j=1:length(D)
    Dr(j)=D(j)/(N*kappas(i)+D(j));
end
Drs{i}=Dr;
sizes(i)=size(Drs{i},1);

%% Try direct lingam
clc
clear all
rng(32)
X = randn(2, 100);
X = [X(1, :); X(1, :)];
A = Dlingam(X); 

%%
X = csvread("../lingamX2.csv",1);
X = X';

tic
A = Dlingam(X);
toc
