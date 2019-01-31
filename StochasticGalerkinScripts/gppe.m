function [ gp ] =gppe( elems,nodelist)
%calculates the gauss point of the triangle in physical domain
%   Detailed explanation goes here
s=size(elems);
gp=zeros(s(1,1),3);
for i=1:s(1,1)
    n1=elems(i,2);
    n2=elems(i,3);
    n3=elems(i,4);
    [nx1 ny1]=find(nodelist(:,1)==n1);
    [nx2 ny2]=find(nodelist(:,1)==n2);
    [nx3 ny3]=find(nodelist(:,1)==n3);
    x1=nodelist(nx1,2);
    y1=nodelist(nx1,3);
    x2=nodelist(nx2,2);
    y2=nodelist(nx2,3);
    x3=nodelist(nx3,2);
    y3=nodelist(nx3,3);
    gp(i,1)=elems(i,1);
    gp(i,2)=(1/3)*(x1+x2+x3);
    gp(i,3)=(1/3)*(y1+y2+y3);
end

end

