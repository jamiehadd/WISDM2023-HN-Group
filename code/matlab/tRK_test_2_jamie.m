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


[~,its] = tRK(W,Y_true,X_0,num_its);

%record errors
errs_prod = [];
for j = 1:num_its+1
    est = its{j} - X_true;
    errs_prod = [errs_prod,norm(est(:))];
end

%plot errors vs iterations
semilogy(0:num_its,errs,'b-',0:num_its,errs_prod,'r--','LineWidth',3)
% hold on
% Create ylabel
ylabel('$\| \mathbf{\mathcal{X}}^{(t)} - \mathbf{\mathcal{X}}^\ddagger \|^2_F$',...
    'FontSize',22,...
    'Interpreter','latex');

% Create xlabel
xlabel('iteration $k$','FontSize',22,'Interpreter','latex');
legend('FacTRK using $\mathbf{\mathcal{U}}, \mathbf{\mathcal{V}}$','TRK using $\mathbf{\mathcal{A}}$','Interpreter','latex')
%title('T-Product Tensor RK')
set(gca,'FontSize',18)

