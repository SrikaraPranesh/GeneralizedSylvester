function [ z ] = GenSylvesterJacobiPrecon( M1,M2,p)
% Application of preconditioner for large matrices
%   Detailed explanation goes here


sm2=size(M2);
sm1=size(M1);

z=zeros(sm1(1,1),sm2(1,2));

p=p*(M2);

for i=1:sm2(1,2)
    p(:,i)=(p(:,i)) ;
    z(:,i)=M1\p(:,i);
end



end

