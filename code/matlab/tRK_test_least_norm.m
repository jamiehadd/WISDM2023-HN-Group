addpath('tproduct toolbox 2.0 (transform)/')

%% Generate Tensors - T-Product
l = 15;
p = 7;
n = 10;
m = 12;

%generate tensors
A = randn(m,l,n);
X_true = randn(l,p,n);

%generate consistent measurements
B = tprod(A,X_true);

%run some iterations of tRK
num_its = 10000;
[~,its] = tRK(A,B,zeros(l,p,n),num_its);

A_t = tran(A);
A_prod_inv = tinv(tprod(A,A_t));
X_ln = tprod(tprod(A_t,A_prod_inv),B);

%record errors
errs = [];
for j = 1:num_its+1
    est = its{j} - X_ln;
    errs = [errs,norm(est(:))];
end

%plot errors vs iterations
semilogy(errs)

