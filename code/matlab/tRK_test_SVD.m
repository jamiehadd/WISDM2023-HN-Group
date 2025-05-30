addpath('tproduct toolbox 2.0 (transform)/')

%% Generate Tensors - T-Product

% Comparing random factorization with QR and SVD
% For SVD, mutliply U by S or S by V to obtain new 2-term factorization  

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
[U1, S1, V1] = tsvd(A,transform,'econ');

%Anew = iprod(U1,V1);

%U1 = tprod(U1, S1);
V1 = tran(V1);
V1 = tprod(S1,V1);
%U1 = tprod(U1, S1);

%Anew = tprod(U1,V1);
%disp(Anew - A)

transform.L = @fft; transform.l = p; transform.inverseL = @ifft;
[Q,R] = tqr(A,transform,'econ');

%W=tprod(U,V);

X_0 = randn(m_1,n,p);

Y_true = tprod(A, X_true);

num_its = 2000;

[~,its_X,its_Z] = tRK_t_prod(U,V,Y_true,X_0,num_its);
[~,its_X1,its_Z1] = tRK_t_prod(U1,V1,Y_true,X_0,num_its);
[~,its_X2,its_Z2] = tRK_t_prod(Q,R,Y_true,X_0,num_its);

errs = [];
errs1 = [];
errs2 = [];
for j = 1:num_its+1
    est = its_X{j} - X_true;
    errs = [errs,norm(est(:))];
    est1 = its_X1{j} - X_true;
    errs1 = [errs1,norm(est1(:))];
    est2 = its_X2{j} - X_true;
    errs2 = [errs2,norm(est2(:))];
end
%close all
%plot errors vs iterations
%size(errs)
%errs
figure;
semilogy(1:num_its+1,errs,'r',1:num_its+1,errs1,'b',1:num_its+1,errs2,'g')
hold on
title('T-Product Tensor RK')
xlabel('iterations')
ylabel('error')
legend('A = UV factorization', 'A = SVD factorization', 'A = QR factorization')

