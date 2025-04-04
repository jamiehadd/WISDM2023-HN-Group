function [H,T] = conv2tensor(sizeI,h)
%% Rewrite 2d convolution as a t-product
% Inputs:
% sizeI: the input image size [m1, n1]
% h: m2 x n2 convolution kernel
% Outputs:
% H: n x n1 x m tensor
% T: (mn) x (m1n1) doubly block-circulant matrix

% Example: find the transformation matrix for conv2
% m1 = 2; n1 = 3; m2 = 2; n2 = 2; 
% m = m1+m2-1; %3
% n = n1+n2-1; %4
% I = randn(m1,n1);
% h = randn(m2,n2);
% x = conv2(I,h);
% [H,T] = conv2tensor(I,h);
% I2 = I(end:-1:1,:)'; % m1 x n1
% y = reshape(T*I2(:),n,m)';  y = y(end:-1:1,:);
% norm(x-y,'fro');

% I2 = [I2(end:-1:1,:); zeros(m-m1,n1)];
% I2ten = zeros(n1,1,m); I2ten(:,1,:) = I2';
% z = tprod(H,I2ten)'; % n x 1 x m
% z = z(end:-1:1,:);
% norm(x-z,'fro')

m1 = sizeI(1);
n1 = sizeI(2);
[m2,n2] = size(h);
m = m1+m2-1;
n = n1+n2-1;
h2 = zeros(m,n);
h2(m-m2+1:m,1:n2) = h;
H = zeros(n,n1,m);


for i = 1:m
    temp = circulant(h2(m-i+1,:)');%,1);
    H(:,:,i) = temp(:,1:n1);    
end


% Generate the doubly block-circulant matrix
Q = cell(m,1);
for i = 1:m
    Q{i} = sparse(H(:,:,i));
end

vc = (1:m)';
vr = [1; (m:-1:(m-m1+2))'];
T = cell2mat(Q(toeplitz(vc,vr)));