addpath('tproduct toolbox 2.0 (transform)/')

%% Generate Tensors - T-Product
m_2 = 10;
m_1 = 5;
p = 30;
n = 7;
m = 50;

X_true = randn(m_1,n,p);
A = randn(m,m_1,p);

X_0 = randn(m_1,n,p);

Y_true = tprod(A, X_true);

num_its = 2000;


[~,its] = tRK(A,Y_true,X_0,num_its);

%record errors
errs_trk = [];
for j = 1:num_its+1
    %est = tprod(A,its{j} - X_true);
    est = its{j} - X_true;
    errs_trk = [errs_trk,norm(est(:))];
end

[~,its_trgs] = tRGS(A,Y_true,X_0,num_its);

%record errors
errs_trgs = [];
for j = 1:num_its+1
    %est = tprod(A,its_trgs{j} - X_true);
    est = its_trgs{j} - X_true;
    errs_trgs = [errs_trgs,norm(est(:))];
end

%plot errors vs iterations
semilogy(0:num_its,errs_trk,'b-',0:num_its,errs_trgs,'r--','LineWidth',3)
% hold on
% Create ylabel
ylabel('$\|\mathbf{\mathcal{A}}(\mathbf{\mathcal{X}}^{(k)} - \mathbf{\mathcal{X}}) \|^2_F$',...
    'FontSize',22,...
    'Interpreter','latex');

% Create xlabel
xlabel('iteration $k$','FontSize',22,'Interpreter','latex');
legend('TRK','TRGS','Interpreter','latex')
%title('T-Product Tensor RK')
set(gca,'FontSize',18)

