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

A = tprod(U,V);
transform = struct;
transform.L = @fft; transform.l = p; transform.inverseL = @ifft;
[U1,V1] = tqr(A,transform,'econ');

%W=tprod(U,V);

X_0 = randn(m_1,n,p);

Y_true = tprod(A, X_true);

num_its = 2000;

[~,its_X,its_Z] = tRK_t_prod(U,V,Y_true,X_0,num_its);
[~,its_X1,its_Z1] = tRK_t_prod(U1,V1,Y_true,X_0,num_its);

errs = [];
errs1 = [];
for j = 1:num_its+1
    est = its_X{j} - X_true;
    errs = [errs,norm(est(:))];
    est1 = its_X1{j} - X_true;
    errs1 = [errs1,norm(est1(:))];
end
close all
%plot errors vs iterations
%size(errs)
%errs
semilogy(1:num_its+1,errs,'r',1:num_its+1,errs1,'b')
hold on
title('T-Product Tensor RK')
xlabel('iterations')
ylabel('error')
legend('A = UV factorization', 'A = QR factorization')

