function [Uf, Lambdaf] = Nystrom_eps_pinv_single(A, n, k, l, rngseed)
%%% Function that produces the Nystrom approximation using the epsilon
%%% pseudoinverse when all the steps are computed in single precision

%%% Parameters %%%
%n : size of matrix A
%k : rank of low rank approx
%l : oversampling parameter
%rngseed: Fixed seed for reproducibility

epsilon = 1.19e-7;

%%% Nystrom method using the epsilon-pseudoinverse %%%
A=single(A);
rng(rngseed);
omega=single(randn(n,k+l));

[Q,~]=qr(omega,0);

Y=A*Q;

B=Q'*Y;
[U,Sigma]=eig((B+B')/2);

%Shift for the epsilon pseudoinverse
epsilon = epsilon*norm(Y,'Fro');
search=diag(Sigma);
idx=find(abs(search)>epsilon);
U=U(:,idx);Sigma=sqrt(Sigma(idx,idx));

F=Y*U/Sigma;
[Uf,Lambdaf,~] = svd(F,0);
Lambdaf=Lambdaf^2;