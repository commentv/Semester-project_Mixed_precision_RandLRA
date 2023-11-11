clear;

%Script that produces the figure that plots the absolute error when the 
%mixed precision nystrom cholesky algorithm is applied on a chosed test 
%problem

example = 'poldecay'; % choose the problem
n = 1e2; % size of A
number_trials=10; %To take the mean of the error
mvp_vec = ['d','s','h']; %Vector for the precisions
k_vec = 1:2:50; %Span of the rank of the approximation
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

Err_matrix_half = zeros(1,length(k_vec));
Err_matrix_simple = zeros(1,length(k_vec));
Err_matrix_double = zeros(1,length(k_vec));
count=0;
    
A = create_example(example,n,q,G,ksi,p);
sing_decay = diag(A);
sing_decay = sing_decay(2:51);

%Perform the approximation in the different precisions%
for j = 1:length(k_vec)
    k = k_vec(j);
    for m = 1:number_trials

        [U,Lambda] = Nystrom(A,n,k,l,'h',m);
        Err_matrix_half(j) = Err_matrix_half(j) + norm(A-U*Lambda*U');
        
        [U,Lambda] = Nystrom(A,n,k,l,'s',m);
        Err_matrix_simple(j) = Err_matrix_simple(j) + norm(A-U*Lambda*U');
    
        [U,Lambda] = Nystrom(A,n,k,l,'d',m);
        Err_matrix_double(j) = Err_matrix_double(j) + norm(A-U*Lambda*U');
    end

    Err_matrix_half(j) = Err_matrix_half(j)/number_trials;
    
    Err_matrix_simple(j) = Err_matrix_simple(j)/number_trials;

    Err_matrix_double(j) = Err_matrix_double(j)/number_trials;

    count=count+1
end
%Plot%
figure;
ax_1 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'log');
title(ax_1,sprintf('%s, %s, cholesky',example,param))
ylabel(ax_1,'$\|A-\hat{A}_{k}\|_{2}$','Interpreter','latex')
xlabel(ax_1,'k');
axis(ax_1,[1e0 50 1e-4 1])
hold(ax_1,'on')
semilogy(ax_1,1:50,sing_decay,'--k');%Plot the best low rank approximation%
semilogy(ax_1,k_vec,Err_matrix_double,'-*',"MarkerEdgeColor",[0.8500 0.3250 0.0980],...
    "MarkerFaceColor",[0.8500 0.3250 0.0980]);
semilogy(ax_1,k_vec,Err_matrix_simple,'-*',"MarkerEdgeColor",[0.9290 0.6940 0.1250],...
    "MarkerFaceColor",[0.9290 0.6940 0.1250]);
semilogy(ax_1,k_vec,Err_matrix_half,'-*',"MarkerEdgeColor",[0.4940 0.1840 0.5560],...
    "MarkerFaceColor",[0.4940 0.1840 0.5560]);
legend(ax_1,'SVD','double','single','half','Location','southwest');