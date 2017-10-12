%Test_ANN_KL

clear
%-----------------------------------------------------------------------

%Assign parameter
new_dim_start=7;    
new_dim_end  =7;
kfold=10;                       %Defining cross-validation
iteration=10;                    %Defining test times

num_train=100;                  %Defining the input
num_test=1;

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

    %---------------------------------------------------------------
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    X_whole=[X;X_star];
    Y_whole=[Y;Y_starorig];

    [X_whole,X_digi]=DataDigiRescale(X_whole);
    X=X_whole(1:np_train,:);
    X_star=X_whole(np_train+1:end,:);

    [Y_whole,Y_digi]=DataDigiRescale(Y_whole);
    Y=Y_whole(1:np_train,:);
    Y_starorig=Y_whole(np_train+1:end,:);

    Y_starorig=ones(np_test,1)*(10.^ Y_digi).*Y_starorig;

    %--------------------------------------------------------------

    X=X';
    Y=Y';
    X_star=X_star';
    Y_starorig=Y_starorig';

    %------------------------------------------------------------------------------------------
    % Testing

    for new_dim=new_dim_start:new_dim_end
        
        %options = struct('ker','poly','arg',[2,1],'new_dim',new_dim); 
        %options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
        options = struct('ker','sigmoid','arg',[0.0000000000000001,0],'new_dim',new_dim); %0.000001
        [Y_star_kpca,t_kpca]=ANN_KPCA(X,Y,X_star,options,kfold); 
        [Y_star_svd,t_svd]=ANN_SVD(X,Y,X_star,new_dim,kfold);
        
        Y_star_kpca=(ones(np_test,1)*(10.^ Y_digi).*Y_star_kpca')';
        Y_star_svd=(ones(np_test,1)*(10.^ Y_digi).*Y_star_svd')';

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
