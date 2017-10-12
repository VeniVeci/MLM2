%Test_PreImage
%Comparing reconstruction performance between different metohds and
%parameters
clear

%Assign parameter
num=160;
dim=10;
new_dim_end=5;

%-----------------------------------------------------------------------
%temp variable
Index_Fig=1;
j=1;

%----------------------------------------------------------------
%Define Input

%Generating test data
%{
    X_temp=rand(1,num);
    for i =1:dim
    X(i,:)=X_temp.^i;
    end
%}
    [X,Y,X_star,Y_starorig]=DataGenerate(num,1);
    X=Y';

%------------------------------------------------------------
%Testing
%options = struct('ker','rbf','arg',5000,'new_dim',new_dim)

for i = 1:5
    options = struct('ker','poly','arg',[2,0],'new_dim',5); 
    %options = struct('ker','rbf','arg',50,'new_dim',3); 
    %options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',3); %0.000001
    Rec_SSE(:,j)=PreImageTest_KPCA(X,options,new_dim_end);
    j=j+1;
    
end


% options = struct('ker','rbf','arg',[50,0],'new_dim',3);
% Rec_SSE(:,1)=PreImageTest_KPCA(X,options,new_dim_end);
% 
% options = struct('ker','poly','arg',[2,0],'new_dim',3);
% Rec_SSE(:,2)=PreImageTest_KPCA(X,options,new_dim_end);
% 
% options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',3);
% Rec_SSE(:,3)=PreImageTest_KPCA(X,options,new_dim_end);


Rec_SSE