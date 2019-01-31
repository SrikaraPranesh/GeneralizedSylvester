function [MV] = GSylvester_MATVEC(K,G,rv_dom,X)
% Computes MATVEC of generalized sylvester eqn operator
%   Detailed explanation goes here
sk=size(K{1,1});
sg=size(G{1,1});

MV=zeros(sk(1,1),sg(1,1));

for i=1:rv_dom
    MV=MV+(K{i,1}*X*G{i,1});
end

end

