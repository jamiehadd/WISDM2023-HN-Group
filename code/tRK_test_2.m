addpath('tproduct toolbox 2.0 (transform)/')

%% Generate Tensors - T-Product
m_2 = 10;
m_1 = 15;
p = 4;
n = 4;
m = 15;

X_true = randn(m_1,n,p);
U = randn(m, m_2,p);
V = randn(m_2,m_1,p);

W=tprod(U,V);

X_0 = zeros(m_1,n,p);

Y_true = tprod(W, X_true);

num_its = 3000;

[~,its_X,~] = tRK_t_prod(U,V,Y_true,X_0,num_its);

close all

errs = [];
lsqerr = [];
for j = 1:num_its+1
    est = its_X{j} - X_true;
    errs = [errs,norm(est(:))];
    est = tprod(W,est);
    lsqerr = [lsqerr, norm(est(:))];
end 

semilogy(lsqerr, 'color','blue') 
hold on

[~,its_X] = tRK(W,Y_true,X_0,num_its);

errs = [];
lsqerr = [];
for j = 1:num_its+1
    est = its_X{j} - X_true;
    errs = [errs,norm(est(:))];
    est = tprod(W,est);
    lsqerr = [lsqerr, norm(est(:))];
end 

semilogy(lsqerr, 'color','red')

title('T-Product Tensor RK')




