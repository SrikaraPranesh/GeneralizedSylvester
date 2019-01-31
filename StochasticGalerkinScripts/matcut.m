function K = matcut(K_uc,D)
% Generates the constrained stiffness matrix
% inputs are unconstrained stiffness matrix and column or row id to be
% removed. DOF number to be removed must be provided in increasing order

s_d=size(D);
for i=1:s_d(1,1)
    
    if (i==1)
        [m,n] = size(K_uc);
        z=D(1,1);
        d1 = K_uc(1:z-1,1:z-1);
        d2 = K_uc(1:z-1,z+1:n);
        d3 = K_uc(z+1:m,1:z-1);
        d4 = K_uc(z+1:m,z+1:n);
        K = [d1 d2; d3 d4];
    else
        D=D-1;
        z1= D(i,1) ;   
        [m,n] = size(K);
        d1 = K(1:z1-1,1:z1-1);
        d2 = K(1:z1-1,z1+1:n);
        d3 = K(z1+1:m,1:z1-1);
        d4 = K(z1+1:m,z1+1:n);
        K = [d1 d2; d3 d4];
    end
end
    
        
end

