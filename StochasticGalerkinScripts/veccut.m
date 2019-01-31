function F = veccut(F_uc,D)
% To remove member i from C and return D with size 1 less.
s_d=size(D);
for i=1:s_d(1,1)
    if (i==1)
       [m,n] = size(F_uc);
       z=D(1,1);
       d1 = F_uc(1:z-1);
       d2 = F_uc(z+1:m);
       F = [d1; d2];
    else
        D=D-1;
        z1= D(i,1) ;
        [m,n] = size(F);
        d1 =F(1:z1-1);
        d2 =F(z1+1:m);
        F = [d1; d2];
    end
end
    
   
   
   
   
   
   
   
   
   
   
end
