%Test_ANN_KLF

clear

%-----------------------------------------------------------------------
%Assign parameter
new_dim_start=3;    
new_dim_end  =3;
kfold=10;                       %Defining cross-validation
iteration=1;                    %Defining test times

num_train=120;                  %Defining the input
num_test=1;

%-----------------------------------------------------------------------
%temp variable
Index_Fig=1;
Rec=zeros(1,9);

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

%------------------------------------------------------------------------------------------
% Testing
    tic;
    Y_star_full= ANN_prediction_crossvalidation(X,Y,X_star,kfold);
    t_full=toc

    for new_dim=new_dim_start:new_dim_end
        
        %options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
        %options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
        options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
        [Y_star_kpca,t_kpca]=ANN_KPCA(X,Y,X_star,options,kfold); 
        [Y_star_svd,t_svd]=ANN_SVD(X,Y,X_star,new_dim,kfold);

        %------------------------------------------------------------
        %Record
        SquErr_kpca=(Y_starorig-Y_star_kpca).^2;
        RateErr_kpca=(abs(Y_starorig-Y_star_kpca))./Y_starorig;
        RecSsErr(new_dim,1)=sum(SquErr_kpca(:));                   %Square sum err
        RecARsErr(new_dim,1)=sum(RateErr_kpca(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(new_dim,1)=t_kpca;                                 %time
    
        SquErr_svd=(Y_star_svd-Y_starorig).^2;
        RateErr_svd=(abs(Y_star_svd-Y_starorig))./Y_starorig;
        RecSsErr(new_dim,2)=sum(SquErr_svd(:));                   %Square sum err
        RecARsErr(new_dim,2)=sum(RateErr_svd(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(new_dim,2)=t_svd;                                 %time
        
        SquErr_full=(Y_starorig-Y_star_full).^2;
        RateErr_full=(abs(Y_starorig-Y_star_full))./Y_starorig;
        RecSsErr(1,3)=sum(SquErr_full(:));                   %Square sum err
        RecARsErr(1,3)=sum(RateErr_full(:))/(np_test*Dim_Y); %Relative sum err
        RecTime(1,3)=t_full;                                 %time


        %Plot
        figure(Index_Pic)
        %title(sprintf('NonRescaled New dim= %g', new_dim))
        plot(Y_star_kpca,'-ob')
        hold on
        plot(Y_star_svd,'-+r')
        plot(Y_star_full','g')
        plot(Y_starorig,'-dk')
        legend('Y^*-KPCA','Y^*--LPCA','Y^*--full','Y^*--real','Location','northeast')
        hold off
        Index_Pic=Index_Pic+1;
        
    end
    
    Record=[RecSsErr,RecARsErr,RecTime];
    Rec=[Rec;Record];

    
end
    
Rec=removerows(Rec,'ind',1);    
