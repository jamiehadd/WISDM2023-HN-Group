function [X_reorder] = reorder_tensor(X, sizeX)
%--------------------------------------------------------------------------
% [X_reorder] = reorder_tensor(X, sizeX)
%
% This routine reorders the tensor to be used in the deblurring
% algorithm--see paper.
%
% Inputs:
% X         -- tensor two be reordered
% sizeX    
%
% Outputs       
% X_reorder -- resized tensor
%--------------------------------------------------------------------------
% Authors: WiSDM 2023 Tensor Group (A. Castillo, J. Haddock, Iryna Hartsock, 
%          P. Hoyos, L. Kassab, A. Kryshchenko, K. Larripa, D. Needell, S.
%          Suryanarayanan, K. Yacoubou Djima)
% 
% Last Modified: May 6, 2025
%--------------------------------------------------------------------------

    % Variable initialization
    m = sizeX(1);
    n = sizeX(2);
    p = sizeX(3);
    X_reorder = zeros(n,p,m);
        X_unfold = zeros(m*n,p);
    for i = 1:p
        X_unfold(1:m*n,i) = reshape(X(:,:,i),[],1);
    end 
    startpt = 1;
    endpt = n;
    
    % Reordering
    for i = 1:m
        X_reorder(1:n,1:p,i) = X_unfold(startpt:endpt,1:p);
        startpt = startpt+n;
        endpt = endpt+n;
    end
end