function MX = mvmx(X)
%% mvmx forms the matrix-vector multiplication M * X where M = A + mu*B + Z*Z'
%  Inputs:
%    X: the projected matrix in Rayleigh Ritz method
%  Outputs:
%    MX: the matrix-vector multiplication
%%
global gA gB gZ gMu
MX = gA * X + gMu * (gB * X) + gZ * (gZ' * X);
end