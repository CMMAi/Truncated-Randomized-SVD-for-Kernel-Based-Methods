function [U,S,V] = rsvd0(X,r)
% Step 1: Sample column space of X with P matrix
ny = size(X,2);
rng(0);
P = randn(ny,r);
Z = X*P;
[Q,R] = qr(Z,0);
% Step 2: Compute SVD on projected Y=Q’*X;
Y = Q'*X;
[UY,S,V] = svd(Y,'econ');
U = Q*UY;
