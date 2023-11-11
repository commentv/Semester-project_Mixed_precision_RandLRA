clear;

%Script that produces the figure when both the both mixed precision nystrom
%algorithms are applied to the RBFK matrix derived from the unifrom
%distribution

%%% Parameters %%%
n=500; %Size of the matrix
k_vec=[1:15 16:2:30]; % rank of low rank approx
l=0; %Oversampling parameter
mvp_vec = ['d','s','h']; %mvp : precision of matrix-matrix multiplication : 'd' (double), 's' (simple), 'h' (half)
rngseed = 1;
sigmaTest=1; %Sigma parameter for matrix abalone

rng(rngseed);
Uniform = rand(1,n);

%Form the RBFK matrix from the uniform vector
Test=zeros(n,n);
for i = 1:n
    for j = 1:n
        Test(i,j)=exp(-abs(Uniform(i)-Uniform(j))^2/sigmaTest^2);
    end
end
[U,Sigma,~] = svd(Test);
Sigma = diag(Sigma);
best_approx = Sigma(2:31);

figure;
ax_1 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'log');
title(ax_1,sprintf('Uniform, \\sigma=%g, n=%i, cholesky',sigmaTest,n))
ylabel(ax_1,'$\|A-\hat{A}_{k}\|_{2}$','Interpreter','latex')
xlabel(ax_1,'k');
axis(ax_1,[1 30 1e-16 1e3])
hold(ax_1,'on')

figure;
ax_2 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'log');
title(ax_2,sprintf('Uniform, \\sigma=%g, n=%i, eps pinv',sigmaTest,n))
ylabel(ax_2,'$\|A-\hat{A}_{k}\|_{2}$','Interpreter','latex')
xlabel(ax_2,'k');
axis(ax_2,[1 30 1e-16 1e3])
hold(ax_2,'on')

semilogy(ax_1,1:30,best_approx,'--k');
semilogy(ax_2,1:30,best_approx,'--k');%Plot the best low rank approximation%

%Perform the approximation for both methods
i=0;
for mvp = mvp_vec
    err_vec_nys_d = [];
    err_vec_nys_pinv_d = [];
    err_vec_nys_d_F = [];
    err_vec_nys_pinv_d_F = [];
    for k = k_vec
            [U,lambda] = Nystrom(Test,n,k,l,mvp,rngseed);

            mat_err = Test-U*lambda*U';

            err=norm(mat_err);
            err_vec_nys_d = [err_vec_nys_d err];

            [U,lambda] = Nystrom_eps_pinv(Test,n,k,l,mvp,rngseed);
    
            mat_err = Test-U*lambda*U';
    
            err=norm(mat_err);
            err_vec_nys_pinv_d = [err_vec_nys_pinv_d err];

          i=i+1
    end
    semilogy(ax_1,k_vec,err_vec_nys_d,'-*');
    semilogy(ax_2,k_vec,err_vec_nys_pinv_d,'-*');
end

legend(ax_1,'SVD','double','single','half','Location','southwest');
legend(ax_2,'SVD','double','single','half','Location','southwest');