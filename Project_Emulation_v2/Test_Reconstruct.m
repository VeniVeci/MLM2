% Test_Reconstruct

clear

%Assign parameter
new_dim_start=4;    
new_dim_step =1;
new_dim_end  =5;

num_train=200;                  %Defining the input parameters
num_test=300;

%Assign memory for recoder
SSE=[];

%------------------------------------------------------------------------------------
%Defining the input

[X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

% Y=Y*10000;
% Y_starorig=Y_starorig*10000;

%------------------------------------------------------------------------------------------
% Testing------------------------------

%-------------------------------------------------------------------------
% SVD Reconstruct
for new_dim=new_dim_start:new_dim_step:new_dim_end

    [Y_star]=SVD_reconstruct(Y,Y_starorig,new_dim);
    SE =(Y_star-Y_starorig).^2;
    SSE_SVD(new_dim,1)=sum(SE(:));

end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% KPCA Reconstruct
Y_starorig=Y_starorig';     %Get dataset ready (in the right format)
Y=Y';
for new_dim=new_dim_start:new_dim_step:new_dim_end
  % options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
    options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
  % ptions = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
    [Y_star]=Kpca_reconstruct(Y,Y_starorig,options,10);
    SE =(Y_star-Y_starorig).^2;
    SSE_KPCA(new_dim,1)=sum(SE(:));

end

%-------------------------------------------------------------------------

% SSE=[SSE,SSE_SVD];

SSE=[SSE_KPCA,SSE_SVD];
    
    
    