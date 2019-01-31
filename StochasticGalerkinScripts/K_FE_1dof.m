function [K] = K_FE_1dof(local_nodes,elems,c)
%global stiffness matrix for CST
%   Detailed explanation goes here

n_ele=size(elems);
n_nodes=size(local_nodes);
K=zeros(n_nodes(1,1),n_nodes(1,1));
for i=1:n_ele(1,1)
    n1=elems(i,2);
    n2=elems(i,3);
    n3=elems(i,4);
    node1=find(local_nodes(:,1)==n1);
    node2=find(local_nodes(:,1)==n2);
    node3=find(local_nodes(:,1)==n3);
    temp=n_nodes(1,1);
    R=zeros(3,temp);
    t1_n1=node1;
    R(1,t1_n1)=1;
    t1_n2=node2;
    R(2,t1_n2)=1;
    t1_n3=(node3);
    R(3,t1_n3)=1;
    nodes=[local_nodes(node1,2) local_nodes(node1,3);local_nodes(node2,2) local_nodes(node2,3);local_nodes(node3,2) local_nodes(node3,3)];
    x1=nodes(1,1);
    y1=nodes(1,2);
    x2=nodes(2,1);
    y2=nodes(2,2);
    x3=nodes(3,1);
    y3=nodes(3,2);
    y23=y2-y3;
    y31=y3-y1;
    y12=y1-y2;
    x32=x3-x2;
    x13=x1-x3;
    x21=x2-x1;
    t=0.01;       % Element thickness
    B=(t/2)*[y23  y31  y12 ; x32  x13  x21];
%     D=(E/(1-((nu)^2)))*[1 nu 0;nu 1 0;0 0 ((1-nu)/2)];
    K_2=c(i,1)*(B')*B;
    K_1=R'*(K_2)*R;
    K=K+K_1;
end
end

