function [X, lam, resHist] = lobpcgsr(A,B,Z,k,tol,mu,maxOut,maxIn)
% lobpcgsr solves A x = lambda B x where A and B are semi-definite
% and the common null-space is Z.
% Inputs:
%   A: constrained graph Laplacian
%   B: constrained graph Laplacian
%   k: LOBPCG block size
%   tol: tolerance of LOBPCG
%   mu: parameter of spectral transformation, used in method (2) and (3)
%   maxOut: maximum number of LOBPCG iterations
%   maxIn: maximum number of inner PCG iterations
% Outputs:
%   X: eigenvectors corresponding to the k smallest eigenvalues
%   lam: k smallest eigenvalues
%   resHist: residual history

global gA gB gZ gMu gMaxIn
gA = A; gB = B; gZ = Z; gMu = mu; gMaxIn = maxIn;
n = size(A,1);
X0 = sprand(n,k,0.1);
t = tic;
[X,sig,~,~,resHist] = lobpcg(X0,-gB,@mvmx,@innerpcg,[],tol,maxOut);
lam = -1 ./ sig - mu;
disp(['Running time:' num2str(toc(t)) ' seconds']);
disp(['Number of iterations:' num2str(size(resHist, 2))]);

