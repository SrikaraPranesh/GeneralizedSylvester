function [F] = f_2d_1dof(elems,nodelist,local_nodes,q)
%force vector
% force vector for uniformly distributed surface traction and its assembly
% for global force vector
   % load intensity 'This can be changes based on requirement'
n_ele=size(elems);
n_node=size(local_nodes);
F=zeros(n_node(1,1),1);

for i=1:n_ele(1,1)
    f_temp=zeros(n_node(1,1),1);
    n1=elems(i,2);
    n2=elems(i,3);
    n3=elems(i,4);
    node1=find(local_nodes(:,1)==n1);
    node2=find(local_nodes(:,1)==n2);
    node3=find(local_nodes(:,1)==n3);   
    t1_n1=node1;
    t1_n2=node2;
    t1_n3=node3;
    X=[nodelist(n1,2);nodelist(n2,2);nodelist(n3,2)];
    Y=[nodelist(n1,3);nodelist(n2,3);nodelist(n3,3)];
    A = polyarea(X,Y);
    f_temp(t1_n1,1)=(q*A)/3;
    f_temp(t1_n2,1)=(q*A)/3;
    f_temp(t1_n3,1)=(q*A)/3;
    F=F+f_temp;
        
end
end

