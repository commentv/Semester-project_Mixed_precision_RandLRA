function [U, Lambda] = Nystrom_single(A, n, k, l, rngseed)
%%% Function that produces the Nystrom approximation using the epsilon
%%% pseudoinverse when all the steps are computed in single precision

%%% Parameters %%%
%n : size of matrix A
%k : rank of low rank approx
%l : oversampling parameter
%rngseed: Fixed seed for reproducibility

epsilon = 1.19e-7;

%%% Nystrom method using the cholesky factorization%%%
A = single(A);
rng(rngseed);
omega=single(randn(n,k+l));

[Q,~]=qr(omega,0);

Y=A*Q;

%Shift to make the cholesky factorisation stable
mu = epsilon*norm(Y,'fro');

Y_mu = Y+mu*Q;
B=Q'*Y_mu;
C=chol((B+B')*0.5);
F=Y_mu/C;
[U,Sigma,~] = svd(F,0);

%Extract the shift
Lambda = max(0,Sigma^2-mu*eye(k+l));