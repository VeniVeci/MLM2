%Test_ANN_KL

clear
%-----------------------------------------------------------------------

%Assign parameter
new_dim_start=1;    
new_dim_end  =5;
kfold=10;                       %Defining cross-validation
iteration=10;                    %Defining test times

num_train=40;                  %Defining the input
num_test=200;

%-----------------------------------------------------------------------
%temp variable
Index_Fig=1;
Rec=zeros(1,6);

%-----------------------------------------------------------------------
%Getting input
for i = 1:iteration

[X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

    %X=load('Inputs2.txt');
    %Y=load('Outputs2.txt');
    %X_star=load('InputTest14');
    %Y_starorig=load('OutputTest14');
    
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

X=X';
Y=Y';
X_star=X_star';
Y_starorig=Y_starorig';

%------------------------------------------
%Make dataset positive defined
% [Dim_Y,Num_Y]=size(Y);
% Y=Y+ones(Dim_Y,Num_Y)*2;
% [Dim_Y_starorig,Num_Y_starorig]=size(Y_starorig);
% Y_starorig=Y_starorig+ones(Dim_Y_starorig,Num_Y_starorig)*2;

%------------------------------------------------------------------------------------------
% Testing

    for new_dim=new_dim_start:new_dim_end
        
        %options = struct('ker','poly','arg',[1.4,10],'new_dim',new_dim); 
        %options = struct('ker','rbf','arg',[10,0],'new_dim',new_dim); %[50,0]
        
%       options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
        options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
        
        [Y_star_kpca,t_kpca]=ANN_KPCA(X,Y,X_star,options,kfold); 
        [Y_star_svd,t_svd]=ANN_SVD(X,Y,X_star,new_dim,kfold);

        %------------------------------------------------------------
        %Record
        SquErr_kpca=(Y_starorig-Y_star_kpca).^2;
        RateErr_kpca=(abs(Y_starorig-Y_star_kpca))./abs(Y_starorig);
        RecSsErr(new_dim,1)=sum(SquErr_kpca(:));                   %Square sum err
        RecARsErr(new_dim,1)=sum(RateErr_kpca(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(new_dim,1)=t_kpca;                                 %time

        SquErr_svd=(Y_star_svd-Y_starorig).^2;
        RateErr_svd=(abs(Y_star_svd-Y_starorig))./abs(Y_starorig);
        RecSsErr(new_dim,2)=sum(SquErr_svd(:));                   %Square sum err
        RecARsErr(new_dim,2)=sum(RateErr_svd(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(new_dim,2)=t_svd;                                 %time

        %Plot
        figure(Index_Fig)
        title(sprintf('NonRescaled New dim= %g', new_dim))
        plot(Y_star_kpca,'-ob')
        hold on
        plot(Y_star_svd,'-+r')
        plot(Y_starorig,'-dk')
        hold off
        legend('Y^*-KPCA.','Y^*--LPCA','Y^*--real.','Location','northeast')
        Index_Fig=Index_Fig+1;
        
    end
    
    Record=[RecSsErr,RecARsErr,RecTime];
    Rec=[Rec;Record];
 
end

Rec=removerows(Rec,'ind',1);    


temp=(sum(Rec'))';
%temp=temp';
index_void=find(temp==0);

Rec2=removerows(Rec,'ind',index_void);    
        