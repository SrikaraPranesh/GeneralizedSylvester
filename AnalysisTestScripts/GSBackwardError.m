function  GSBackwardError(siv,cnv,ns,flag,testnum)
%GSBACKWARDERROR compares the backward error of 2 term generalized Sylvester
%equation with relative residual.
%   siv -- An integer indicating the size of component matrices.
%   cnv -- A 1X4 vactor of integers indicating the condition number of
%           the component matrices.
%   ns  -- Number of seeds to be used to generate the matrices
%   flag -- use flag=1 to generate the plots, else 0

sv = [1:1:ns]';
dv = zeros(length(sv),1);
msz = 500; %size of the markers in scatter plot
msz1 = 30; %size of the markers in the legend
for i=1:length(sv)
    rng(sv(i,1))
    n = siv;
    
    A1=gallery('randsvd',n,cnv(1,1));
    B2=gallery('randsvd',n,cnv(1,4));
    
    A2=gallery('randsvd',n,cnv(1,3));
    B1=gallery('randsvd',n,cnv(1,2));
    
    A=(kron(B1',A1))+kron(B2',A2);
    b=randn((n*n),1);
    F=reshape(b,n,n);
    
    x=A\b;
    X=reshape(x,n,n);
    cnsol(i,1) = cond(X);
    
    res=norm(b-(A*x),inf);
    r=b-(A*x);
    
    nbe(i,1)=(res)/(norm(b,inf)+(norm(A,inf)*norm(x,inf)));
    
    T1=[(norm(A1,'fro')*kron((B1'*X'),eye(n))) (norm(A2,'fro')*kron((B2'*X'),eye(n)))];
    T2=[(norm(B1,'fro')*kron(eye(n),(A1*X)))  (norm(B2,'fro')*kron(eye(n),(A2*X)))];
    T3=norm(F,'fro')*eye((n*n));
    H=[T1 T2 -T3];
    
    pH = pinv(H);
    abe(i,1) = norm((pH*r));
    
    %%%%% Estimating the condition number
    Cnumber(i,1) = sqrt(5)*norm((A\H))/norm(X,'fro');
    
    numer = 2*((norm(A1,'fro')*norm(B1,'fro'))+(norm(A2,'fro')*norm(B2,'fro')))*norm(X,'fro')...
        + norm(F,'fro');
    Cnumber(i,2) = sqrt(5)*norm((A\eye(length(A))))*(numer/norm(X,'fro'));
    
    fprintf('sample number %d out of %d\n',i,length(sv));
    
end


%%%%%%% generate plots %%%%
if (flag == 1)
    
    %%%% Scatter plots of backward errors
    scatter(sv,abe,msz,'d');
    hold on
    scatter(sv,nbe,msz,'o');
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',40)
    xlabel('Seed')
    [~, objh] = legend('actual backward error','relative residual');
    objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
    set(objhl, 'Markersize', msz1); %// set marker size as desired
    
    
    str_e1 = sprintf('%0.1e',cnv(1,1));
    str_e2 = sprintf('%0.1e',cnv(1,2));
    str_e3 = sprintf('%0.1e',cnv(1,3));
    str_e4 = sprintf('%0.1e',cnv(1,4));
    title(sprintf('Case %d',testnum));
    
    %%%% Scatter plot of condition number of the solution matrix
    scatter(sv,cnsol,msz,'d');
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',40)
    xlabel('Seed')
    ylabel('condition number')
    title(sprintf('Case %d',testnum));
    
    
    %%%% Scatter plot of the condition number
    scatter(sv((1:30),1),Cnumber((1:30),1),msz,'d');
    hold on
    scatter(sv((1:30),1),Cnumber((1:30),2),msz,'o');
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',40)
    xlabel('Seed')
    ylabel('condition number')
    title(sprintf('Case %d',testnum));
    
    [~, objh] = legend('strong condition number','weak condition number');
    objhl = findobj(objh, 'type', 'patch'); %// objects of legend of type line
    set(objhl, 'Markersize', msz1); %// set marker size as desired
    
end

end

