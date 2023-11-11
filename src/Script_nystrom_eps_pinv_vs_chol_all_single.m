clear;

%Script that produces the figure that compares the absolute error when both
%single precision nystrom algorithms are applied to a chosed test problem

example = 'expdecay'; % choose the problem
n = 1e2; % size of A
number_trials=10; %To take the mean of the error
k_vec = 1:1:30;%Span of the rank of the approximation
l=0; %oversampling parameter

%%% Parameters for the exponential decay %%%
q = 1; % rate of exp. decay: 0.1 slow, 0.25 med, 1 fast

%%% Parameters for the psd noise %%%
G = randn(n);% To get same G and ksi on the measures
ksi = 1e-1; % 1e-4; 1e-2, 1e-1

%%% Parameters for the polynomial decay %%%
p = 2; % rate of decay: 0.5 slow, 1 med, 2 fast

%%% Setting parameter strings for the title of the plot %%%
switch example
    case 'expdecay'
            param = sprintf('q=%g, n=%i',q,n);

       case 'psdNoise'
            param = sprintf('ksi=%g, n=%i',ksi,n);

       case 'poldecay'
            param = sprintf('p=%g, n=%i',p,n);

    case 'stairdecay'
            param = sprintf('n=%i',n);
end

Err_matrix_simple = zeros(1,length(k_vec));
Err_matrix_all_simple = zeros(1,length(k_vec));
count=0;

A = create_example(example,n,q,G,ksi,p);
sing_decay = diag(A);
sing_decay = sing_decay(2:31);

%Perform the approximation when everything is in single precision for the
%method using the cholesky factorization and the espilon pseudoinverse
for j = 1:length(k_vec)
    k = k_vec(j);
    for m = 1:number_trials
        
        [U,lambda] = Nystrom_single(A,n,k,l,m);
        Err_matrix_simple(j) = Err_matrix_simple(j) + norm(A-U*lambda*U');
        [U,lambda] = Nystrom_eps_pinv_single(A,n,k,l,m);
        Err_matrix_all_simple(j) = Err_matrix_all_simple(j) + norm(A-U*lambda*U');
    
    end
    Err_matrix_simple(j) = Err_matrix_simple(j)/number_trials;
    Err_matrix_all_simple(j) = Err_matrix_all_simple(j)/number_trials;

    count=count+1
end
%Plot%
figure;
ax_1 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'log');
title(ax_1,sprintf('%s, %s, cholesky vs eps pinv',example,param))
ylabel(ax_1,'$\|A-\hat{A}_{k}\|_{2}$','Interpreter','latex')
xlabel(ax_1,'k');
axis(ax_1,[1e0 30 1e-8 1e0])
hold(ax_1,'on')
semilogy(ax_1,1:30,sing_decay,'--k');%Plot the best low rank approximation%
semilogy(ax_1,k_vec,Err_matrix_simple,'-or');
semilogy(ax_1,k_vec,Err_matrix_all_simple,'-xb');
legend(ax_1,'SVD','cholesky all single','eps pinv all single','Location','northeast');

