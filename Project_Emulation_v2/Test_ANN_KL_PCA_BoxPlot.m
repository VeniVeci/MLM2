%Test_ANN_KL-PCA_BoxPlot

clear

%Assign parameter
new_dim_start=1;    
new_dim_step =1;
new_dim_end  =5;

num_train=80;                  %Defining the input parameters
num_test=100;

kfold=10;                       %Defining cross-validation

%------------------------------------------------------------------------------------
%Defining the input

[X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

X=X';
Y=Y';
X_star=X_star';

% Y=Y*10000;
% Y_starorig=Y_starorig*10000;

%------------------------------------------------------------------------------------------
% Testing------------------------------

    for new_dim=new_dim_start:new_dim_step:new_dim_end
        
        %options = struct('ker','poly','arg',[2,0],'new_dim',new_dim); 
        options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
%         %options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
%         options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);

        [Y_star_kpca,t_kpca]=ANN_KPCA(X,Y,X_star,options,kfold); 
        [Y_star_svd,t_svd]=ANN_SVD(X,Y,X_star,new_dim,kfold);

        Y_star_kpca=Y_star_kpca';
        Y_star_svd=Y_star_svd';
        
        %------------------------------------------------------------ 
        %Record
             
%         Y_star_kpca=real(Y_star_kpca);
        
        %KPCA part
        SquErr=(Y_starorig-Y_star_kpca).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_kpca)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_kpca(:,new_dim)   =SSErr;
        RecRateSsErr_kpca(:,new_dim)=RatSsErr;
        RecRateErr_kpca(:,new_dim)=RateErr;
        RecTime_kpca(:,new_dim)=t_kpca;
  
        %SVD part
        SquErr=(Y_starorig-Y_star_svd).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_svd)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_svd(:,new_dim)   =SSErr;
        RecRateSsErr_svd(:,new_dim)=RatSsErr;
        RecRateErr_svd(:,new_dim)=RateErr;
        RecTime_svd(:,new_dim)=t_svd;
                     
    end
    
%     save('RecSSErr_kpca_nptrain100.txt','RecSSErr_kpca','-ascii')
%     save('RecSSErr_svd_nptrain100.txt','RecSSErr_svd','-ascii')
    
figure(1)
boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of KPCA-ANN'));

figure(2)
boxplot(RecRateSsErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Square Sum Error Rate of Each pixel of LPCA-ANN'));

figure(3)
boxplot(RecRateErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Average Error Rate of Each pixel of KPCA-ANN'));

figure(4)
boxplot(RecRateErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
title(sprintf('Average Error Rate of Each pixel of LPCA-ANN'));


% 
% 
%         Mass=reshape(Y_star_svd,100,100);
%         Mass=Mass';
%         figure21=figure('InvertHardcopy','off','Color',[1 1 1]);
%         axes1 = axes('Parent',figure21,'FontSize',28,'FontName','Times');        
%         box(axes1,'on');
%         hold(axes1,'all');
%         contourf(x2,x1,Mass,'Parent',axes1)
%         % Create xlabel
%         xlabel('x / m','FontSize',28,'FontName','Times');
%         % Create ylabel
%         ylabel('y / m','FontSize',28,'FontName','Times');
%         title(sprintf('(b)'),'FontSize',28,'FontName','Times')
%         colorbar('peer',axes1,'FontSize',28,'FontName','Times')


dim_physical=1;
nptx=100;
Lx=10;
Ly=10;

dx=Lx/(nptx-1);
dy=Ly/(nptx-1);
x1=0:dx:Lx;
x2=0:dy:Ly;

n=43;

Mass=reshape(Y_starorig(n,:),nptx,nptx);
figure(9)

contourf(x2,x1,Mass)
        % Create xlabel
        xlabel('\xi / cm','FontSize',24,'FontName','Times');
        % Create ylabel
        ylabel('\eta / cm','FontSize',24,'FontName','Times');
        title(sprintf('(a)'),'FontSize',24,'FontName','Times')
        colorbar('FontSize',24,'FontName','Times')
% contourf(Mass)
% title(sprintf('Real'));

Mass=reshape(Y_star_kpca(n,:),nptx,nptx);
figure(10)
contourf(x2,x1,Mass)
        % Create xlabel
        xlabel('\xi / cm','FontSize',24,'FontName','Times');
        % Create ylabel
        ylabel('\eta / cm','FontSize',24,'FontName','Times');
        title(sprintf('(b kpca)'),'FontSize',24,'FontName','Times')
        colorbar('FontSize',24,'FontName','Times')

Mass=reshape(Y_star_svd(n,:),nptx,nptx);
figure(11)
contourf(x2,x1,Mass)
        % Create xlabel
        xlabel('\xi / cm','FontSize',24,'FontName','Times');
        % Create ylabel
        ylabel('\eta / cm','FontSize',24,'FontName','Times');
        title(sprintf('(c lpca)'),'FontSize',24,'FontName','Times')
        colorbar('FontSize',24,'FontName','Times')


