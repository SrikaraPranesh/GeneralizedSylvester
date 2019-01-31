% This script reproduces Table 2 of 
% 'Solution of the Sylvester Matrix Equation 
% AXB^T+CXD^T=E -- Judith D Gardiner, Alan J Laub,
% James J Amato, and Cleve B Moler, Gardiner, 
% ACM Transactions on Mathematical Software (TOMS) 
% 18.2 (1992): 223-231. 

clear all
close all

m = 10; n = 4; p = [0;10;20;30;40];
Um = tril(ones(m));
Un = tril(ones(n));

for i=1:length(p)

A = diag([1:1:m]')+Um;
B = eye(n)+(((2^(-p(i,1))))*Un');
C = eye(m)+(((2^(-p(i,1))))*Um');
D = ((2^(-p(i,1)))*eye(n)) - diag([n:-1:1]') + Un;

G = kron(B,A)+kron(D,C);
xact = ones(length(G),1);
Xact = ones(m,n);
b = double(mp(mp(G,34)*mp(xact,34),34));
E = reshape(b,m,n);

x = G\b;
X = reshape(x,m,n);

%%% Forward error and relative residual
NE = norm(X-Xact,inf)/norm(Xact,inf);

R = ((A*X*B')+(C*X*D')-E);
NR = norm(R,'fro')/(((norm(X,'fro'))*((norm(A,'fro')*norm(B,'fro'))...
     +(norm(C,'fro')*norm(D,'fro'))))+norm(E,'fro'));

%%% Actual backward error
T1 = [(norm(A,'fro')*kron((B*X'),eye(m))) (norm(C,'fro')*kron((D*X'),eye(m)))];
T2 = [(norm(B','fro')*kron(eye(n),(A*X))) (norm(D','fro')*kron(eye(n),(C*X)))];
T3 = norm(E,'fro')*eye((m*n));
H = [T1 T2 -T3];

r = reshape(R,(m*n),1);
abe = norm(((H)\r));

result(i,:) = [NR abe];

Cnumber(i,1) = sqrt(3)*norm((G\H))/norm(X,'fro');
Cnumber(i,2) = cond(G);

fprintf('Sample number %d\n',i)

end


