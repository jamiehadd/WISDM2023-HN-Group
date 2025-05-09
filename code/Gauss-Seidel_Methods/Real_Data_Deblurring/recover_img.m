function [Z_rec] = recover_img(Z,sizeX)
%--------------------------------------------------------------------------
% [Z_rec] = recover_img(Z,sizeX)
%
% This routine recovers a tensor for the deblurring algorithm--see paper.
%
% Inputs:
% Z         -- tensor two be recovered
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
    Z_rec = zeros(m,n,p);
    Z_unfold = zeros(n*m,p);
    startpt = 1;
    endpt = n;

    % Recovering
    for i = 1:m
        Z_unfold(startpt:endpt,1:p) = Z(:,:,i);
        startpt = startpt+n;
        endpt = endpt+n;
    end
    
    for i = 1:p
        Z_rec(:,:,i) = reshape(Z_unfold(:,i),m,n);
    end
end
