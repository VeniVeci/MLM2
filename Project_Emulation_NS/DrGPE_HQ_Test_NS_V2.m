% DrGPE_HQ_Test_NS_V2
% A HQ terminal for dimension reduction GPE _Special design for Navier
% stock equation for the special format of dataset( e.g. using U,V,P as Y
% which are multi-output.
% In V2 all result are saved
%% Initialize
% clear
% close all 

% load('Data5kPCAinfExact.mat')
preoptions.type='Exp';  %'Exp' method is more stable and non-parametric
%  preoptions.neighbor=10;

num_test=100;
Test_split=201;


mean_SE=[];
MSSE_U=[];
MSSE_V=[];
MSSE_P=[];

SSE=[];
SSE_U=[];
SSE_V=[];
SSE_P=[];

for i=1:length(num_train_ref)
    
   
    % Special design for special structure of Diffusion maps.
    if strcmp(DrMethod,'DiffusionMaps' )
        dim_extract=dim_new+1;
    else 
        dim_extract=dim_new;
    end
    
    % Collect dataset
    num_train=num_train_ref(i);
    Z=Z_Rec(1:num_train,1:dim_extract,i);
    model=model_Rec(i);
    
    Y=Y_Data(1:num_train,:);
    X=X_Data(1:num_train,:);
    
    X_star=X_Data(Test_split:Test_split+num_test-1,:);
    Y_starorig=Y_Data(Test_split:Test_split+num_test-1,:);
    
    % Solving Pre-image with increasing dimensions
    
    [num_Z,dim_Z]=size(Z);    
    %Clean memory
    m=zeros(num_test,dim_Z);
    s=zeros(num_test,dim_Z);
      
    for j=1:dim_Z
            hyp=hyp_Rec(i,j);
            [m(:,j) s(:,j)] = gp(hyp, InfMethod, meanfunc, covfunc, likfunc, X, Z(:,j), X_star);  
    end
        
    %Clean memory
    SSE=zeros(num_test,dim_new);
    
    for k=1:dim_new
        
        switch DrMethod    
            case 'kPCA'
                Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);
            case 'DiffusionMaps'
                Y_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
            case 'Isomaps'

            otherwise 
                error('No such DR method')
        end

        Rec_Y_star(:,:,k,i)=Y_star;

    end
    
end



for i=1:length(num_train_ref)
    for k=1:dim_new
        
        Y_star=Rec_Y_star(:,:,k,i);
        
        U_star=Y_star(:,1:Dim_U);
        V_star=Y_star(:,1+Dim_U:Dim_U+Dim_V);
        P_star=Y_star(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
        
        U_starorig=Y_starorig(:,1:Dim_U);
        V_starorig=Y_starorig(:,1+Dim_U:Dim_U+Dim_V);
        P_starorig=Y_starorig(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
             
        
        SquErr=(Y_starorig-Y_star).^2;
        SquEr_U=(U_starorig-U_star).^2;
        SquEr_V=(V_starorig-V_star).^2;
        SquEr_P=(P_starorig-P_star).^2;

        
        SSE(:,k,i)=sum(SquErr');
        SSE_U(:,k,i)=sum(SquEr_U');
        SSE_V(:,k,i)=sum(SquEr_V');
        SSE_P(:,k,i)=sum(SquEr_P');
    end
end
        
        
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

for i=1:length(num_train_ref)

     VarRelaError=abs(s./m);
     sumVarRelaError(i,:)=sum(VarRelaError);
%     figure(10+i)
%     boxplot(VarRelaError)
%     title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
%     set(gca,'yscale','log');
%     xlabel('Index of principal direction')
%     ylabel('Predicted variance ratio')
               
    figure(i)
    boxplot(SSE(:,:,i))
    title(sprintf('Boxplot of Square sum error for different subspace dimension -%s -%d Training points', DrMethod, num_train_ref(i)))
    set(gca,'yscale','log');
    xlabel('Dimensions')
    ylabel('SSE')

%      mean_SE(i,:)=mean(SSE);     
%      MSSE_U(i,:)=mean(SSE_U);
%      MSSE_V(i,:)=mean(SSE_V);
%      MSSE_P(i,:)=mean(SSE_P);
     
end













%BOXPLOT

markers = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
[num_temp,dim_temp]=size(mean_SE);
figure
hold on
for i=1:num_temp
    semilogy(mean_SE(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('Mean of Square sum error for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('MSSE')

%----------------------------
[num_temp,dim_temp]=size(MSSE_U);
figure
hold on
for i=1:num_temp
    semilogy(MSSE_U(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('MSSE of U field for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('MSSE')

%----------------------------
[num_temp,dim_temp]=size(MSSE_V);
figure
hold on
for i=1:num_temp
    semilogy(MSSE_V(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('MSSE of V field for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('MSSE')

%----------------------------
[num_temp,dim_temp]=size(MSSE_P);
figure
hold on
for i=1:num_temp
    semilogy(MSSE_P(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('MSSE of P field for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('MSSE')

%ACUTAL PLOT
% Dim=5;
% yplot=SSE(:,Dim);
% h=boxplot(yplot);
% obj = findobj(h,'tag','Outliers');
% ydata=get(obj,'ydata');
% 
% Value_median=median(yplot);
% 
% [yplot_sort,index]=sort(yplot);






% index=[17:20];
% figure
% hold on
% plot(Y_starorig(index,:)','b-')
% plot(Y_star(index,:)','k--')
% hold off


    
