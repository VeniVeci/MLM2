% Test prediction on PCA-GPR and KPCA-GPR. 
%
% Description:
% This is a modificated version that do comparation and show the best
% result.
%
% steps for the program:
% 1.run PCA-GPR and KPCA-GPR on given dataset for desinated principal
% components. pick up best result and record them.
% 2. run it many times and packup the result.
%
% Library needed
% 1.Gaussian Process
% 2.DimensionReduction
% 3.Data_Demonstration
%
% Modifications:
% 4-Sep-2013, WeiX, first edition 
%---------------------------------------------------------------------

%--------------------------------------
% Add path (Similar to including library in C)
path_current=pwd;
path_GP = strcat(path_current,'/Gaussian Process');
path_DR = strcat(path_current,'/DimensionReduction');
addpath(path_GP,path_DR);


clear

%--------------------------------------
%Assign parameter
new_dim_start=1;    
new_dim_end  =7;

num_train=40;                  %Defining the input
num_test=1;

iteration=100;                   %Defining test times

%-----------------------------------------------------------------------
%temp variable
Rec=zeros(1,6);

%-----------------------------------------------------------------------

for i = 1:iteration  
    
    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

        %X=load('Inputs2.txt');
        %Y=load('Outputs2.txt');
        %X_star=load('InputTest14');
        %Y_starorig=load('OutputTest14');

    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);
    
    %Rescale 
    Y=Y*0.001;
    Y_starorig=Y_starorig*0.001;
    
    %------------------------------------------------------------------------------------------
    % Testing------------------------------

    for new_dim=new_dim_start:new_dim_end
        
        options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
        %options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
        %options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
        
        [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
        [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);

        %------------------------------------------------------------
        %Record
        SquErr_kpca=(Y_starorig-Y_star_kpca).^2;
        RateErr_kpca=(abs(Y_starorig-Y_star_kpca))./Y_starorig;
        RecSsErr(new_dim,1)=sum(SquErr_kpca(:));
        RecARsErr(new_dim,1)=sum(RateErr_kpca(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(new_dim,1)=t_kpca;

        SquErr_svd=(Y_star_svd-Y_starorig).^2;
        RateErr_svd=(abs(Y_star_svd-Y_starorig))./Y_starorig;
        RecSsErr(new_dim,2)=sum(SquErr_svd(:));
        RecARsErr(new_dim,2)=sum(RateErr_svd(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(new_dim,2)=t_svd;

%         %Plot
%         figure(new_dim)
%         title(sprintf('NonRescaled New dim= %g', new_dim))
%         plot(Y_star_kpca,'-ob')
%         hold on
%         plot(Y_star_svd,'-+r')
%         plot(Y_starorig,'-dk')
%         hold off
%         legend('Y^*-KPCA.','Y^*--LPCA','Y^*--real.','Location','northeast')
        
    end
    
    Record=[RecSsErr,RecARsErr,RecTime];
    
    %-------------------------------------
    % Searching for best result. recording with principal components and time
    
    RecSsErr=real(RecSsErr);    %Not sure if it is right to do this here. Real number are supposed to be here but it is not always the case, so this line is added to ensure real number.
    
    RecSsErr(~RecSsErr)=NaN;    %ensure the minimal search will not search for result of zero value.

    [SsErr_min_value_kpca,SsErr_min_PCs_kpca]=min(RecSsErr(:,1));
    SsErr_min_time_kpca=RecTime(SsErr_min_PCs_kpca,1);
    
    [SsErr_min_value_pca,SsErr_min_PCs_pca]=min(RecSsErr(:,2));
    SsErr_min_time_pca=RecTime(SsErr_min_PCs_pca,2);
    
    Record_Kpca=[SsErr_min_PCs_kpca,SsErr_min_value_kpca,SsErr_min_time_kpca];   %formate=[used PCs, SSE, Time]
    Record_pca =[SsErr_min_PCs_pca, SsErr_min_value_pca, SsErr_min_time_pca];    %formate=[used PCs, SSE, Time]
    Record=[Record_Kpca,Record_pca];
    
    %--------------------------
    %Recording
    Rec=[Rec;Record];
    
end
    
Rec=removerows(Rec,'ind',1);    

%-----------------------------------------------
% Deleting blank line. It help showing a brief result
 temp=(sum(Rec'))';  
 temp=temp';
 index_void=find(temp==0);

 Rec2=removerows(Rec,'ind',index_void);
 
 %------------------------------------------------

 % Plot the Time comparation without boundary
 bound_down =0;
 bound_up   =100000000;
 Array_HistCompare(Rec2(:,3),Rec2(:,6),bound_down,bound_up);   
 
 % Plot the SSE comparation with boundary
 bound_down =0;
 bound_up   =100;
 Array_HistCompare(Rec2(:,2),Rec2(:,5),bound_down,bound_up);     

