clear;

%Script that produces the figure when both the both mixed precision nystrom
%algorithms are applied to the RBFK matrix derived from the abalone data

%Extract the abalone data
A=readtable('..\data\abalone.txt');
Abalone=table2array(A(:,2:end));

%%% Parameters %%%
k_vec=[1 10:10:150]; % rank of low rank approx
l=0; %Oversampling parameter
mvp_vec = ['d','s','h']; %mvp : precision of matrix-matrix multiplication : 'd' (double), 's' (simple), 'h' (half)
rngseed = 1;

nAb = size(Abalone,1); %Number of rows of Abalone
dAb = size(Abalone,2); %Number of columns of Abalone
sigmaAb=1; %Sigma parameter for matrix abalone

%Form the RBFK matrix from the abalone data
AbaloneD=zeros(nAb,nAb);
for i = 1:nAb
    for j = 1:nAb
        norm_diff = vecnorm(Abalone(i,:)-Abalone(j,:));
        AbaloneD(i,j)=exp(-(norm_diff^2)/sigmaAb^2);
    end
end

[U,Sigma,~] = svd(AbaloneD);
Sigma = diag(Sigma);

best_approx_d = Sigma(2:151);

figure;
ax_1 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'log');
title(ax_1,sprintf('Abalone Dense, \\sigma=%g, cholesky',sigmaAb))
ylabel(ax_1,'$\|A-\hat{A}_{k}\|_{2}$','Interpreter','latex')
xlabel(ax_1,'k');
axis(ax_1,[1 150 1e-1 1e3])
hold(ax_1,'on')

figure;
ax_3 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'log');
title(ax_3,sprintf('Abalone Dense, \\sigma=%g, eps pinv',sigmaAb))
ylabel(ax_3,'$\|A-\hat{A}_{k}\|_{2}$','Interpreter','latex')
xlabel(ax_3,'k');
axis(ax_3,[1 150 1e-1 1e3])
hold(ax_3,'on')

semilogy(ax_1,1:150,best_approx_d,'--k');
semilogy(ax_3,1:150,best_approx_d,'--k');%Plot the best low rank approximation%

%Perform the approximations for both methods
i=0;
for mvp = mvp_vec
    err_vec_nys_d = [];
    err_vec_nys_pinv_d = [];
    for k = k_vec
        [U,lambda] = Nystrom(AbaloneD,nAb,k,l,mvp,rngseed);
        mat_err = AbaloneD-U*lambda*U';
        err=norm(mat_err);
        err_vec_nys_d = [err_vec_nys_d err];

        [U,lambda] = Nystrom_eps_pinv(AbaloneD,nAb,k,l,mvp,rngseed);
        mat_err = AbaloneD-U*lambda*U';
        err=norm(mat_err);
        err_vec_nys_pinv_d = [err_vec_nys_pinv_d err];

          i=i+1
    end
    semilogy(ax_1,k_vec,err_vec_nys_d,'-*');
    semilogy(ax_3,k_vec,err_vec_nys_pinv_d,'-*');
end

legend(ax_1,'SVD','double','single','half','Location','southwest');
legend(ax_3,'SVD','double','single','half','Location','southwest');