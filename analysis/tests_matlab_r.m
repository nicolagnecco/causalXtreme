%% Test mentappr
clear all;

x = csvread("temp_csv/x.csv", 1);
mentappr(x)

%% Test pwling
clear all;

X1 = csvread("temp_csv/X1.csv", 1);
pwling(X1, 1)


%% Test direct lingam
clear all;

X3 = csvread("temp_csv/X3.csv", 1)';

[~, ~, ~, K] = pwDling(X3)

X4 = csvread("temp_csv/X4.csv", 1)';

[~, ~, ~, K] = pwDling(X4)


%% Test pwlingc
clc;
X5 = csvread("temp_csv/X5.csv", 1)';
tic
[~, ~, ~, K] = pwDling(X5)
toc
 