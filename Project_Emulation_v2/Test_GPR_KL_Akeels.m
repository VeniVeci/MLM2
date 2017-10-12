%Test_GPR_KL

clear

%Assign parameter
new_dim_start=1;    
new_dim_step =2;
new_dim_end  =9;

dim_physical=1;
nptx=100;
Lx=10;
Ly=10;

dx=Lx/(nptx-1);
dy=Ly/(nptx-1);
x1=0:dx:Lx;
x2=0:dy:Ly;
num_train=40;                  %Defining the input parameters
num_test=400;

%------------------------------------------------------------------------------------
%Defining the input

[X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);
%------------------------------------------------------------------------------------------
% Testing------------------------------

    for new_dim=new_dim_start:new_dim_step:new_dim_end
        
        %options = struct('ker','poly','arg',[2,10],'new_dim',new_dim); 
        options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
        %options = struct('ker','rbf','arg',50,'new_dim',new_dim); 
        %options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);
        
        [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
        [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);

        %------------------------------------------------------------ 
        %Record
             
        Y_star_kpca=real(Y_star_kpca);
        
        %KPCA part
        SquErr=(Y_starorig-Y_star_kpca).^2;
        means=abs(mean(Y_starorig,2));
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
        means=abs(mean(Y_starorig,2));
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

% Create figure
figure(3)
% figure3 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure3,'XTickLabel',{'2','4','6','8','10','12'},...
%     'XTick',[1 2 3 4 5 6],...
%     'FontSize',20,...
%     'FontName','Times');
% 
% %set(gca,'XTickLabel',{' '})
boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
% Create xlabel
xlabel('Number of principal components','Units','points','FontSize',20,...
    'FontName','Times');

title(sprintf('(b) KPCA'));

% Create ylabel
ylabel('Relative square error','FontSize',20,'FontName','Times');
figure(4)
% Create figure
% figure4 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure4,'XTickLabel',{'2','4','6','8','10','12'},...
%     'XTick',[1 2 3 4 5 6],...
%     'FontSize',20,...
%     'FontName','Times');

boxplot(RecRateSsErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});

% Create xlabel
xlabel('Number of principal components','Units','points','FontSize',20,...
    'FontName','Times');

title(sprintf('(a) LPCA'));

% Create ylabel
ylabel('Relative square error','FontSize',20,'FontName','Times');

if dim_physical>1
    
n=143;


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
        
        
n=322;


Mass=reshape(Y_starorig(n,:),nptx,nptx);
figure(12)

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
figure(13)
contourf(x2,x1,Mass)
        % Create xlabel
        xlabel('\xi / cm','FontSize',24,'FontName','Times');
        % Create ylabel
        ylabel('\eta / cm','FontSize',24,'FontName','Times');
        title(sprintf('(b kpca)'),'FontSize',24,'FontName','Times')
        colorbar('FontSize',24,'FontName','Times')

Mass=reshape(Y_star_svd(n,:),nptx,nptx);
figure(14)
contourf(x2,x1,Mass)
        % Create xlabel
        xlabel('\xi / cm','FontSize',24,'FontName','Times');
        % Create ylabel
        ylabel('\eta / cm','FontSize',24,'FontName','Times');
        title(sprintf('(c lpca)'),'FontSize',24,'FontName','Times')
        colorbar('FontSize',24,'FontName','Times')
else
    timev=0:14/60:7000/60;
 
    kplot(1)=189;
    kplot(2)=96;
    kplot(3)=227;
    kplot(4)=333;
    kplot(5)=67;

figure(9)
for n=1:5
plot(timev,Y_starorig(kplot(n),:)*10)
hold on
plot(timev,Y_star_svd(kplot(n),:)*10,'r')
end
xlabel('Time / min','FontSize',24,'FontName','Times');
ylabel('Temperature / C','FontSize',24,'FontName','Times');

figure(10)
for n=1:5
plot(timev,Y_starorig(kplot(n),:)*10)
hold on
plot(timev,Y_star_kpca(kplot(n),:)*10,'r')
end
xlabel('Time / min','FontSize',24,'FontName','Times');
ylabel('Temperature / C','FontSize',24,'FontName','Times');

end