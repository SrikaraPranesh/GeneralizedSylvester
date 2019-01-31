%%%% This script generates the plots given in 
%%%% 'Backward error and condition number of generalized
%%%%  Syvester equation -- Srikara Pranesh'

clear all
close all

% Size of component matrices
n = 4;

% Vector of condition number of component matrices
cnv(1,:) = [1e15 1e15 1e15 1e15];
cnv(2,:) = [1e15 1e15 1e2 1e15];
cnv(3,:) = [1e15 1e15 1e2 1e2];
ns = 1000;
flag = 1;


for i=1:3
   CouExCnumber{i,1} = GSBackwardError(n,cnv(i,:),ns,flag);
end

for i = 1:length(CouExCnumber)
    dv{i,1} = ((CouExCnumber{i,1}(:,2)-CouExCnumber{i,1}(:,3))./CouExCnumber{i,1}(:,2))*100;
end



%%%%%%% PRINT THE LATEX TABLE INTO A FILE %%%%%%%

% creating a text file to print the GMRES iteration table
fid1 = fopen('CounterExampleCond.txt','w');
fprintf(fid1,'Column 1 is seed used to generate the matrix \n');
fprintf(fid1,'Columns 2 is the actual condition number \n');
fprintf(fid1,'Column 3 is the condition number corresponding to linear systems \n');
fprintf(fid1,'\n'); fprintf(fid1,'\n');

for i = 1:length(CouExCnumber)
    [row,col] = size(CouExCnumber{i,1});
    for j = 1:row
    t1  =  CouExCnumber{i,1}(j,1); t2 = CouExCnumber{i,1}(j,2);
    t3  =  CouExCnumber{i,1}(j,3);
    fprintf(fid1,'%d & %6.2e & %6.2e \\\\ \n',t1,t2,t3);
    end
    fprintf(fid1,'\\hline \n');
end

fclose(fid1);


