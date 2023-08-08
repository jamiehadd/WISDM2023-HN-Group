addpath('tproduct toolbox 2.0 (transform)/')

%% Generate Tensors - T-Product
m_2 = 10;
m_1 = 5;
p = 7;
n = 5;
m = 50;

X_true = randn(m_1,n,p);
U = randn(m, m_2,p);
V = randn(m_2,m_1,p);

W=tprod(U,V);

X_0 = randn(m_1,n,p);

Y_true = tprod(W, X_true);

num_its = 2000;

[~,its_X,its_Z] = tRK_t_prod(U,V,Y_true,X_0,num_its);

errs = [];
for j = 1:num_its+1
    est = its_X{j} - X_true;
    errs = [errs,norm(est(:))];
end
close all
%plot errors vs iterations
semilogy(errs)
hold on
title('T-Product Tensor RK')

