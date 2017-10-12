% Test_GPR_LPCA
%
% This script is written to track down the fail prombel when PCA_GPR is
% applied to the data set: "5050PoleFlow"
% 20-Nov-2013 

clear

%Assign parameter

num_train=200;                  %Defining the input parameters
num_test=100;

%Assign memory for recoder
SSE=[];

%----------------------
% parameters
new_dim=2;

for i = 1:new_dim
    RecZstar{i}=[];
    RecZ{i}=[];
end

X_orig=load('Inputs_5050PloeFlow.txt');
Y_orig=load('Outputs_5050PloeFlow.txt');


    splith=401;
    X_star=X_orig(splith:splith+num_test-1,:);
    Y_starorig=Y_orig(splith:splith+num_test-1,:);
    
i=1;
num_train=51;
ind_train=num_train;
RecFilter=[];   %Record point delete form data set
X_orig_copy=X_orig;


%% PCA-GPR test
for k = 1:1:250
    
    %------------------------------------------------------------------------------------
    %Defining the input

%     [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

    TrianStar=1;

    X=X_orig(TrianStar:TrianStar+num_train-1,:);
    Y=Y_orig(TrianStar:TrianStar+num_train-1,:);
    
%     [num_X_orig,~]=size(X_orig);
%     splith=num_X_orig-num_test+1;
% %     splith=401;
%     X_star=X_orig(splith:splith+num_test-1,:);
%     Y_starorig=Y_orig(splith:splith+num_test-1,:);
% 
%     [np_train,Dim_X]=size(X);
%     [np_test, Dim_Y]=size(Y_starorig);

    %------------------------------------------------------------------------------------------
    % Testing------------------------------
    
    [Z,Key] = DataReduc_SVD(Y,new_dim);
    [Z_star,~]= GPR_prediction(X,Z,X_star);
    
    %-------------------------------------
%     % Recording the process
%     for j=1:new_dim
%          Z_temp=nan*ones(num_train_end,1);
%          Z_temp(1:num_train,1)=Z(:,j);
%          RecZ{j}=[RecZ{j},Z_temp];             %RecZ{j} means Z for jth principal component
%     end
%       
%     for j=1:new_dim                
%          Zstar_temp=Z_star(:,j); % Assign memory
%          RecZstar{j}=[RecZstar{j},Zstar_temp];  
%     end
    %-------------------------------------
    
    Y_star=Z_star*Key;
    
    SE =(Y_star-Y_starorig).^2;
    sse_temp=sum(SE(:));
   
    
    if sse_temp> 500
        X_orig=removerows(X_orig,'ind',num_train);    
        Y_orig=removerows(Y_orig,'ind',num_train);  
        
        X_orig_copy(ind_train,:)=zeros(1,3);       
        RecFilter=[RecFilter;ind_train]; 
    else 
        SSE(i,1)=sse_temp;
        i=i+1;
        num_train=num_train+1;
    end
    
    ind_train=ind_train+1;
    
end


% SSE=[SSE,SSE_SVD];

% SSE=[SSE_KPCA,SSE_SVD];
    
    
    