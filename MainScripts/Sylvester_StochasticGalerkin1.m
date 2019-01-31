% This script performs the following tasks. The Generalised Sylvester
% equation is denoted as GSE.
% 1. Computes the actual backward error of GSE.
% 2. Evaluates the estimate for backward error of GSE.
% 3. Estimates the actual normwise condition number of GSE.
% 4. Estimates the normwise condition number of GSE as a stand
%       linear system of equation.
%
% Karhunen-L\"oeve expansion is performed using the algorithm of
%       Srikara Pranesh and Debraj Ghosh. "Faster computation of the
%       Karhunen–Loève expansion using its domain independence property."
%       Computer Methods in Applied Mechanics and Engineering 285 (2015):
%       125-145.

clear all
close all

% Finite element mesh file
fid=fopen('main.msh');

% Converts FE mesh file to MATLAB readable format
[elems,nodelist] = mesh_node(fid,0); % change to 1 for FE mesh plot

% amplitude of uniformly applied force
q = 1;

% boundary condition applied at maximum and minimum
% x- co-ordinate values
b1 = min(nodelist(:,2));
b2 = max(nodelist(:,2));
b = [b1 b2];

%%%% properties of the random field
mean_k = 200;
var_k = 1000;

% Karhunen-L\"oeve expansion performed using
n1 = 35;n2 = 35;
lc = 2.5;
% %%% performing KL expansion on each domain
X1 = max(nodelist(:,2));
X2 = min(nodelist(:,2));
Y1 = max(nodelist(:,3));
Y2 = min(nodelist(:,3));

gp_d = gppe(elems,nodelist);
[evec_gp,eigvals] = KL_DI(X1,Y1,X2,Y2,n1,n2,gp_d,mean_k,var_k,lc);

%%%%
cnts = [2 2;3 2;4 2;2 3;3 3];
rvc = {'H','P'};
save init_data
%%% Generation of polynomial chaos data %%%%
for jj = 1:2
    if (jj==2)
        clearvars -except rqtys cnumbers pp1
        jj = 2;
    end
    load init_data
    for ii = 1:5
        chaos_order = cnts(ii,2);
        norand = cnts(ii,1);
        rv_dom = norand+1;
        I = multiindex(norand,chaos_order);
        sys = rvc{1,jj};
        M = gpcbasis_triples( {sys, I}, {sys, I}, {sys, I} );
        
        %%%% construction of 'B' matrix
        for i = 1:(norand+1)
            G{i,1} = sparse(M(:,:,i));
        end
        pp = length(G{1,1});
        pp1{jj,1}(ii,1) = pp;
        %%%%
        
        %%%% Identifying the boundary nodes %%%
        sb=max(size(b));
        if (sb==1)
            [D1,remain1,F1] = bc_lst_1dof(nodelist,b(1,1),nodelist);
            D2=[];
            remain2=[];
            F2=[];
        elseif (sb==2)
            [D1,remain1,F1] = bc_lst_1dof(nodelist,b(1,1),nodelist);
            [D2,remain2,F2] = bc_lst_1dof(nodelist,b(1,2),nodelist);
        end
        DD  = union(D1,D2);
        FF=union(F1,F2);
        sd=size(DD);
        
        
        if (sd(1,1)>0)
            sr1=max(size(remain1));
            sr2=max(size(remain2));
            if (sr1~=0 && sr2~=0 )
                remain=intersect(remain1,remain2);
            else
                remain=union(remain1,remain2);
            end
        else
            remain=nodelist;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%% Stiffness matrices
        K_m1 = KL_mode_K_1dof(nodelist,nodelist,elems,eigvals,evec_gp);
        
        
        
        %%%% Force vector %%%
        F_loc = f_2d_1dof(elems,nodelist,nodelist,q);
        
        
        K = cell((norand+1),1);
        for i = 1:(norand+1)
            K{i,1} = matcut(K_m1{i,1},DD);
        end
        
        clear K_m1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F2=veccut(F_loc,DD);
        sr=length(F2);
        F=zeros(sr,pp);
        F(:,1)=F2;
        for i=2:pp
            F(:,i)=zeros(sr,1);
        end
        
        
        %%%%% PCG iteration for final SSFEM system of equation %%%%
        sk=size(K{1,1});
        sg=size(G{1,1});
        %%% Setting up the block Jacobi preconditioner
        M1=(K{1,1}); M2=inv(G{1,1});
        X=zeros(sk(1,1),sg(1,1));
        r=F-(GSylvester_MATVEC(K,G,rv_dom,X));
        z=GenSylvesterJacobiPrecon(M1,M2,r);
        p=z;
        rsold=Mat_Frobenius_ip(r,r);
        f_norm=rsold;
        %%%
        
        for i=1:(sk(1,1)*sg(1,1))          % I will have to check what is the actual number I must give here
            Ap=GSylvester_MATVEC(K,G,rv_dom,p);
            denom=Mat_Frobenius_ip(r,z);
            alpha=denom/Mat_Frobenius_ip(Ap,p);
            X=X+alpha*p;
            r=r-alpha*Ap;
            rsnew=Mat_Frobenius_ip(r,r);
            %%%%%
            cutoff_lim=sqrt(rsnew)/sqrt(f_norm);
            %%%%%
            
            
            if (cutoff_lim)<1e-6
                break;
            end
            
            %%%%%
            z=GenSylvesterJacobiPrecon(M1,M2,r);
            numera=Mat_Frobenius_ip(r,z);
            beta=(numera)/(denom);
            %%%%%
            
            p=z+(beta)*p;
            res_norm(i,1)=rsold;
            rsold=rsnew;
            %     fprintf('iteration number %d || residual norm %f \n',i,cutoff_lim);
        end
            
        %%%%computing actual backward error%
        r=F-(GSylvester_MATVEC(K,G,rv_dom,X));
        I=eye(length(K{1,1}));
        l1=length(K{1,1});
        l2=length(G{1,1});
        r=reshape(r,(l1*l2),1);
        H=[];
        denom = 0;
        for i=1:rv_dom
            T=(norm(K{i,1},'fro')*kron(eye(l2,l1),(X*G{i,1})))+(norm(G{i,1},'fro')*kron((X'*K{i,1}),eye(l1,l2)));
            H=[H T];
            denom = denom + (2*(norm(K{i,1},'fro')*norm(G{i,1},'fro')));
        end
        denom = (denom*norm(X,'fro')) + norm(F,'fro');
        
        T=norm(F,'fro')*eye((l1*l2));
        H=[H T];
        pH = pinv(H);
        rs = (norm(r)/denom);
        denom1 = sqrt((((norm(K{1,1},'fro')/norm(inv(G{1,1}),'fro'))*min(svd(X)))^2)+...
            norm(F,'fro')^2);
        mu = denom/denom1;
        
        
        %%%% Estimating various condition numbers
        P = zeros((sr*pp));
        for i = 1:rv_dom
            P = P+kron(G{i,1}',K{i,1});
        end
        Psi = sqrt((2*rv_dom)+1)*(norm(P\H)/norm(X,'fro'));
        nip = min(eig(P));
        if (nip < 0)
           error('matrix is not symmetric and positive definite'); 
        end
        Phi = sqrt((2*rv_dom)+1)*(1/nip)*(denom/norm(X,'fro'));
        lconum = cond(P) +(((1/nip)*norm(F,'fro'))/norm(X,'fro'));
        
        denom11 = (norm(P)*norm(X,'fro')) + norm(F,'fro');
        rs1 = (norm(r)/denom11); 
        rqtys{jj,1}(ii,:) = [norm(pH*r) (mu*rs) rs1];
        cnumbers{jj,1}(ii,:) = [Psi Phi lconum];
        
        fprintf('In %d case out of 5 cases \n',ii);
    end

end
delete init_data.mat;
%%%%%%% PRINT THE LATEX TABLE INTO A FILE %%%%%%%
fid1 = fopen('Result1.txt','w');
fprintf(fid1,'backward errors \n');
fprintf(fid1,'Columns 3, 4 and 5 are for Hermite Polynomials \n');
fprintf(fid1,'Columns 6, 7 and 8 are for Legendre Polynomials \n');
fprintf(fid1,'\n'); fprintf(fid1,'\n');
for i = 1:5
    t0 = pp1{1,1}(i,1);
    t1  =  cnts(i,1); t2 = cnts(i,2);
    t11  = rqtys{1,1}(i,1); t22 = rqtys{1,1}(i,2);
    t3  =  rqtys{1,1}(i,3); t4 = rqtys{2,1}(i,1);
    t33  =  rqtys{2,1}(i,2); t44 = rqtys{2,1}(i,3);
    fprintf(fid1,' & %d &(%d,%d) & %6.2e & %6.2e & %6.2e & %6.2e & %6.2e & %6.2e\\\\ \n',...
        t0,t1,t2,t11,t22,t3,t4,t33,t44);
end

fprintf(fid1,'\n');
fprintf(fid1,'\n');
fprintf(fid1,'condition number \n');
fprintf(fid1,'Columns 3, 4 and 5 are for Hermite Polynomials \n');
fprintf(fid1,'Columns 6, 7 and 8 are for Legendre Polynomials \n');
for i = 1:5
    t0 = pp1{1,1}(i,1);
    t1  =  cnts(i,1); t2 = cnts(i,2);
    t11  = cnumbers{1,1}(i,1); t22 = cnumbers{1,1}(i,2);
    t3  =  cnumbers{1,1}(i,3); t4 = cnumbers{2,1}(i,1);
    t33  =  cnumbers{2,1}(i,2); t44 = cnumbers{2,1}(i,3);
    fprintf(fid1,' & %d &(%d,%d) & %6.2e & %6.2e & %6.2e & %6.2e & %6.2e & %6.2e\\\\ \n',...
        t0,t1,t2,t11,t22,t3,t4,t33,t44);
end

fclose(fid1);
