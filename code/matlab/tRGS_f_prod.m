function [X,its_X,its_Z] = tRGS_f_prod(U,V,Y,X_0,T)
    %record number of row slices
    %m = size(U,1);
    
    [m,l,n] = size(U);
    [m_2,l_2,n_2] = size(V);
    %initialize iterate
    X = X_0;
    Z = tprod(V,X_0);
    its_X = {X};
    its_Z = {Z};
    %iterate
    for t = 1:T
        % %sample row slice from U to update Z in terms of Y residuals
        % i_t = randsample(m,1);
        % U_slice = U(i_t,:,:);
        % Y_slice = Y(i_t,:,:);
        % 
        % %calculate necessary transformations of U_slice
        % U_slice_t = tran(U_slice);
        % U_prod_inv = tinv(tprod(U_slice,U_slice_t));
        % resid_y = tprod(U_slice,Z) - Y_slice;
        % 
        % %RK step
        % Z = Z - tprod(tprod(U_slice_t,U_prod_inv),resid_y);
        % its_Z{end+1} = Z;
%___________________________________________________________________________

        %sample column slice from U to update Z in terms of Y resiuals
        i_t = randsample(l,1);
        U_slice = U(:,i_t,:);
        I = eye(l,l);
        E = zeros(l,1,n);
        E(:,:,1)=I(:,i_t);

        %calculate necessary transformations of U_slice
        U_slice_t = tran(U_slice);
        U_prod_inv = tinv(tprod(U_slice_t,U_slice));
        resid = Y - tprod(U,Z);

        %RGS step
        Z = Z + tprod(tprod(E,tprod(U_prod_inv,U_slice_t)),resid);
        its_Z{end+1} = Z;

%___________________________________________________________________________

        % %sample row slice from V to update X in terms of Z residuals
        % i_t = randsample(m_2,1);
        % V_slice = V(i_t,:,:);
        % Z_slice = Z(i_t,:,:);
        % 
        % %calculate necessary transformations of V_slice
        % V_slice_t = tran(V_slice);
        % V_prod_inv = tinv(tprod(V_slice,V_slice_t));
        % resid_z = tprod(V_slice,X) - Z_slice;
        % 
        % %RK step
        % X = X - tprod(tprod(V_slice_t,V_prod_inv),resid_z);
        % its_X{end+1} = X;
%___________________________________________________________________________

        %sample column slice from V to update X in terms of Z resiuals
        i_t = randsample(l_2,1);
        V_slice = V(:,i_t,:);
        I = eye(l_2,l_2);
        E = zeros(l_2,1,n);
        E(:,:,1)=I(:,i_t);

        %calculate necessary transformations of U_slice
        V_slice_t = tran(V_slice);
        V_prod_inv = tinv(tprod(V_slice_t,V_slice));
        resid = Z - tprod(V,X);

        %RGS step
        X = X + tprod(tprod(E,tprod(V_prod_inv,V_slice_t)),resid);
        its_X{end+1} = X;

    end
end