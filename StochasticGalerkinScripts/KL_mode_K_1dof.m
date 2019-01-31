function [K_sfem] = KL_mode_K_1dof(local_nodes,nodelist,elems,e_value,ev_gp)
%global stiffness matrix for CST with variation in material property
% stiffness corresponding to each KL mode is same through out hence this
% function has been written
% pass the whole eigen function and eigen value set into the function
%   Detailed explanation goes here
% c=*e9;


seval=size(e_value);
for j=1:seval(1,1)
    K_sfem{j,1}=Stiffness(local_nodes,nodelist,elems,e_value,ev_gp,j);
end

end

