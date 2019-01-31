function [D,remain,F] = bc_lst_1dof(nodelist,boundary,local_nodes)
%makes a list of the nodes on the boundary and rows and columns to be cut
%for LST
%   apply BCs for stiffness matrix of FE model with CST element. Boundary
%   must be along a line parellal to y axis. change value of boundaty for
%   various BC's. Can also be used to find interface nodes in FETI if the
%   interface is along a straight line parellal to y-axis.
n_node=size(local_nodes);
a=1;
for i=1:n_node(1,1)
    temp=local_nodes(i,1);
    if(nodelist(temp,2)==boundary)
        B(a,1)=find(local_nodes(i,1)==local_nodes(:,1));  % local node number
        F(a,1)=local_nodes(i,1); % global node
        a=a+1;
    end
end

z=exist('B');
if (z >0)
    D=sort(B);
    %Identifying the remaining nodes and DOF
    remain=setdiff(local_nodes(:,1),F(:,1)); % remaining global nodes
else
    D=[];
    F=[];
    remain=[];
end

end

