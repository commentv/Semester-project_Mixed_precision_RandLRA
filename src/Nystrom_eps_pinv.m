function [Uf, Lambdaf] = Nystrom_eps_pinv(A, n, k, l, mvp, rngseed)
%%% Function that produces the Nystrom approximation using the epsilon
%%% pseudoinverse

%%% Parameters %%%
%n : size of matrix A
%k : rank of low rank approx
%l : oversampling parameter
%mvp : precision of matrix-matrix multiplication : 'd' (double), 's'
%(simple), 'h' (half), 'b' (bfloat16)
%rngseed: Fixed seed for reproducibility

% path to chop (for mixed presicion)
addpath '.\chop-master'

if mvp == 'd'
    epsilon = eps;
elseif mvp == 's'
    epsilon = 1.19e-7;
elseif mvp == 'h'
    epsilon = 9.76e-4;
elseif mvp == 'b'
    epsilon = 3.91e-3;
end

%%% Nystrom method using the epsilon-pseudoinverse %%%
rng(rngseed);
omega=randn(n,k+l);

[Q,~]=qr(omega,0);

%Matrix-matrix product AQ in different precisions
%For detail in the implementation, see the simulation in low precision
%article from Higham
if mvp == 'd'
    Y=A*Q;
elseif mvp == 's'
    A=single(A);
    Q_sing=single(Q);
    Y=single(A*Q_sing);
    Y=double(Y);
else
    opt.format = mvp;
    chop([],opt)
    Y = zeros(n,k+l);
    for i = 1:n
        Y = chop(Y + chop(A(:,i)*Q(i,:)));
    end
end


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