addpath('../tproduct toolbox 2.0 (transform)/')

%% Generate Tensors - T-Product
clear;
clc;

%fixed:
p = 7;
n = 5;

% X_true is 20 x 10 x 30
%m determines whether A = UV is underdetermined or overdetermined
%m_1 determinex whether V is underdetermined or overdetermined
%m_2 is like k in Ma & Needell's paper
m_2 = 10;

%V underdetermined (m_1 > m_2), A overdetermined (m >= 20)--Works
%m = 20;
%m = 30;
%m = 40;
%m_1 = 5;

%V underdetermined (m_1 > m_2), A underdetermined (m < 20)--Doesn't Works
%m = 10;
%m = 15;%
%m_1 = 5;

%V overdetermined (m_1 < m_2), A overdetermined (m > 20)--Doesn't work
% m = 20;
% m = 30;
m = 40;
m_1 = 12;


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