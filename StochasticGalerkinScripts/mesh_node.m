function [elems,nodelist] = mesh_node(fid,fla)
%generates mesh and nodelist file from gmsh .msh file
%   two examples (paste them in command window);
%



C1 = textscan(fid, '%s','delimiter', '\n');
C=C1{1,1}; 
fclose(fid);
s_f=size(C);
nn=str2num(C{5,1});
nln=5+nn;
temp=nln+3;
ne=str2num(C{temp,1});
%eln
% creating nodelist file
fid1=fopen('nodelist.txt','w');
for i=6:nln
    fprintf(fid1,'%s\n',C{i});
end
fclose(fid1);
% creating elems file
fid2=fopen('elems1.txt','w');
for i=(temp+1):(ne+temp)
    temp2=str2num(C{i,1});
    t3=size(temp2);
    if (t3(1,2)==8)
    fprintf(fid2,'%s\n',C{i});
    end
end
fclose(fid2);
load('elems1.txt');
nodelist=load('nodelist.txt');
elems2(:,1)=elems1(:,1);
elems2(:,2:4)=elems1(:,6:8);
fid3=fopen('elems.txt','w');
se=size(elems2);
for i=1:se(1,1)
    t4= num2str(elems2(i,:));
    A{i,1}=t4;
end

for i=1:se(1,1)
    fprintf(fid3,'%s\n',A{i});
end
fclose(fid3);
elems=load('elems.txt');

if (fla==1)
    se=size(elems);
    for i=1:se(1,1)
        n1=elems(i,2);
        n2=elems(i,3);
        n3=elems(i,4);
        l=find(nodelist(:,1)==n1);
        m=find(nodelist(:,1)==n2);
        n=find(nodelist(:,1)==n3);
        x1=nodelist(l,2);
        y1=nodelist(l,3);
        x2=nodelist(m,2);
        y2=nodelist(m,3);
        x3=nodelist(n,2);
        y3=nodelist(n,3);
        X=[x1;x2;x3];
        Y=[y1;y2;y3];
        D=delaunay(X,Y);
        triplot(D,X,Y,'b')
        hold on
    end
end
delete elems1.txt

end

