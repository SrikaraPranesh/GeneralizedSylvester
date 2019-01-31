function [ev_gp,temp_ev] = KL_DI(X1,Y1,X2,Y2,n1,n2,gp,mean_val,var_val,lc)
% THIS IS TOTALLY CORRECTED
% CALCULATES KL CO-EFFICIENTS USING METHOD PROPOSED IN 
% FASTER COMPUTATION OF THE KARHUNEN–LOÈVE EXPANSION USING ITS DOMAIN
% INDEPENDENCE PROPERTY S.PRANESH AND D.GHOSH
% KL_DI(X1,Y1,n1,n2,gp,nm)

% applicable only when FE uses CST elements
%  Refer to the paper for more details
% X1,Y1 are length and width of 2D specimen
% n1 and n2 are number of elements in each direction
% it will have to interpolate eigenfunctions to gauss points
delta_x=abs(X2-X1)/(n1);
delta_y=abs(Y2-Y1)/(n2);
u1=X2+(delta_x/2);
u2=X1-(delta_x/2);
v1=Y2+(delta_y/2);
v2=Y1-(delta_y/2);
cent_x=[u1:delta_x:u2]';
cent_y=[v1:delta_y:v2]';
[cent_posx1,cent_posy1]=meshgrid(cent_x,cent_y);
s1=size(cent_posx1);
s2=size(cent_posy1);
cpx=reshape(cent_posx1,(s1(1,1)*s1(1,1)),1);
cpy=reshape(cent_posy1,(s2(1,1)*s2(1,1)),1);
cp=[cpx cpy];
A=(delta_x*delta_y);
%%%%%% GENERATION OF C and D MATRIX %%%%%%%%%
for i=1:((s1(1,1)*s1(1,1)))
    for j=1:((s1(1,1)*s1(1,1)))
        C(i,j)=var_val*(exp(-((abs(cp(i,1)-cp(j,1))))/(lc)).*(exp(-(abs(cp(i,2)-cp(j,2)))/(lc))))*A;
    end
end
sgp=size(gp);
for i=1:((s1(1,1)*s1(1,1)))
    for j=1:sgp(1,1)
        C1(i,j)=var_val*(exp(-((abs(cp(i,1)-gp(j,2))))/(lc)).*(exp(-(abs(cp(i,2)-gp(j,3)))/(lc))))*A;
    end
end
sc=size(C);
[eigvec1,lambda]=eig(C);
temp_ev2=diag(lambda);
eigvec2=eigvec1;
temp_ev11=temp_ev2;
for i=1:sc(1,1)
    eigvec1(:,i)=eigvec2(:,(sc(1,1)-i+1));
    temp_ev2(i,1)=temp_ev11((sc(1,1)-i+1),:);
end
    

% identifyng dominant KL modes
imp_ratio=0.01;
c=1;
ev_max=max(temp_ev2);
while (temp_ev2(c,1) >= (imp_ratio*ev_max))
    temp_ev1(c,1)=temp_ev2(c,1);
    c=c+1;
end
temp_ev=[1;temp_ev1];

%
s_evec=size(temp_ev1);
for i=1:s_evec(1,1)
    eigvec(:,i)=(eigvec1(:,i)/norm(eigvec1(:,i)));
end

% interpolation of eigenfunctions to gauss points 
ev_gp1=(C1'*eigvec);
for i=1:s_evec(1,2)
    ev_gp1(:,i)=ev_gp1(:,i)/(sqrt(A)*temp_ev1(i,1));
end
sevgp=size(ev_gp1);

tm=ones(sevgp(1,1),1);
ev_gp=[tm*mean_val (ev_gp1)];
% norand=s_evec(1,1);
% nosamp=50000;
% xzi1 = randn(nosamp,norand);
% % 
% xzi=[ones(nosamp,1) xzi1];
% 
% 
% 
% 
% % generation of eigen function plot uncomment only if i want to simulate
% % random process not required for SSFEM computation
% 
% 
% % Z=zeros(s1(1,1),s1(1,1),s_evec(1,2));
% % for i=1:s_evec(1,2)
% %     Z(:,:,i)=reshape(eigvec(:,i),s1(1,1),s1(1,1));
% % end
% % 
% random_process=0;
% % % 
% for i=1:(s_evec(1,1)+1)
%     random_process=random_process+(sqrt(temp_ev(i,1))*ev_gp(1,i)*xzi(:,i));
% end
% % 
% ksdensity(random_process);
% keyboard
% 
% 
% random_process=(random_process/sqrt(A))+200;
% [pdf_mp(:,1) pdf_mp(:,2)]=ksdensity(random_process);








end

