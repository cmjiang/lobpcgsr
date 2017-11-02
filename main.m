clear variables; close all;
addpath('./data');
addpath('./src');

load dataset1 % load matrices A and B
n = size(A,1);

Z = ones(n,1)/sqrt(n); % common null-space

k = 2; % Number of desired eigenvalues

tol = 1e-4; % tolerance

mu = 1e-3; % shift factor in spectral regularization

maxOut = 1000; % max number of lobpcg iterations

maxIn = 4; % max number of inner pcg iterations

% solver eigenproblem A x = lambda B x
[X, lam, resHist] = lobpcgsr(A,B,Z,k,tol,mu,maxOut,maxIn);

% plot the residual history
figure;
semilogy(resHist');