function X=TDMAsolver(A,b)
%% TMDA mthod 
% TDMA method for solving linear equation system

m=length(b);                 % m is the number of rows
X=zeros(m,1);
A(1,2)= A(1,2)  ./ A(1,1);    % Division by zero risk.
b(1)=  b(1)    ./ A(1,1);    % Division by zero would imply a singular matrix

for i=2:m-1
    temp=  A(i,i) - A(i,i-1) .* A(i-1,i);
    A(i,i+1)=  A(i,i+1)  ./ temp;
    b(i)= ( b(i) - A(i,i-1) .* b(i-1) )  ./ temp;
end 

i=m;
X(m)=(b(i) - A(i,i-1) .* b(i-1))  ./ (A(i,i) - A(i,i-1) .* A(i-1,i));

for i=m-1:-1:1
X(i)=  -A(i,i+1) .* X(i+1) + b(i);
end
end