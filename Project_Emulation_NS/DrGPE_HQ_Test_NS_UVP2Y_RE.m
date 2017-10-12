% DrGPE_HQ_Test_NS_UVP2Y_RE
% A HQ terminal for dimension reduction GPE _Special design for Navier
% stock equation for the special format of dataset( e.g. using U,V,P as Y
% which are multi-output. Show relative error.
%% Initialize
% clear
close all 

% load('exp1_kPCA.mat')
preoptions.type='Exp';  %'Exp' method is more stable and non-parametric
% preoptions.neighbor=10;



num_test=300;
Test_split=201;

mean_SE=[];
MRMSSE_U=[];
MRMSSE_V=[];
MRMSSE_P=[];

SSE=[];
SSE_U=[];
SSE_V=[];
SSE_P=[];

for i=3:length(num_train_ref)
    
   
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
    
    for k=5:dim_new
        
        switch DrMethod    
            case 'kPCA'
                Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);
            case 'DiffusionMaps'
                model.options.Ztype=1;
                Y_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
            case 'Isomaps'

            otherwise 
                error('No such DR method')
        end

        
        %Seperate field           
        U_star=Y_star(:,1:Dim_U);
        V_star=Y_star(:,1+Dim_U:Dim_U+Dim_V);
        P_star=Y_star(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
        
        U_starorig=Y_starorig(:,1:Dim_U);
        V_starorig=Y_starorig(:,1+Dim_U:Dim_U+Dim_V);
        P_starorig=Y_starorig(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
             
        %Error record
        SquErr=(Y_starorig-Y_star).^2;
        
        SquEr_U=(U_starorig-U_star).^2;
        SquEr_V=(V_starorig-V_star).^2;
        SquEr_P=(P_starorig-P_star).^2;
        
        MSS_Uorig=mean(U_starorig.^2,2);    %Mean Square Sum of U original field
        MSS_Vorig=mean(V_starorig.^2,2);
        MSS_Porig=mean(P_starorig.^2,2);
        
%         mean_Uorig=mean(U_starorig,2);
%         mean_Vorig=mean(V_starorig,2);
%         mean_Porig=mean(P_starorig,2);
        
%         RSSE_U(:,k)=sum(SquEr_U,2)./MSS_Uorig;          %Relativeu Square Sum Error
%         RSSE_V(:,k)=sum(SquEr_V,2)./MSS_Vorig;
%         RSSE_P(:,k)=sum(SquEr_P,2)./MSS_Porig;        
%         FSRSSE(:,k)=RSSE_U(:,k)+RSSE_V(:,k)+RSSE_P(:,k); %Field Sum Relative Square Sum Error
        
        
        RMSSE_U(:,k)=mean(SquEr_U,2)./MSS_Uorig;        %Relative Mean Square Sum Error
        RMSSE_V(:,k)=mean(SquEr_V,2)./MSS_Vorig;
        RMSSE_P(:,k)=mean(SquEr_P,2)./MSS_Porig;        
        SMRSSE(:,k)=(RMSSE_U(:,k)+RMSSE_V(:,k)+RMSSE_P(:,k))./3; %Field Sum Relative Mean Square Sum Error
        
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

    end
    
     VarRelaError=abs(s./m);
     sumVarRelaError(i,:)=sum(VarRelaError);
%     figure(10+i)
%     boxplot(VarRelaError)
%     title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
%     set(gca,'yscale','log');
%     xlabel('Index of principal direction')
%     ylabel('Predicted variance ratio')
               
    figure(i)
    boxplot(SMRSSE)
    title(sprintf('Boxplot of FSMRSSE -%s -%d Training points', DrMethod, num_train_ref(i)))
    set(gca,'yscale','log');
    xlabel('Dimensions')
    ylabel('SSE')
   
     MRMSSE_U(i,:)=mean(RMSSE_U);    %Mean Sum Relative Square Sum Error
     MRMSSE_V(i,:)=mean(RMSSE_V);
     MRMSSE_P(i,:)=mean(RMSSE_P);
     MMRSSE(i,:)=mean(SMRSSE);  
     
end

%BOXPLOT

markers = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
[num_temp,dim_temp]=size(MMRSSE);
figure
hold on
for i=1:num_temp
    semilogy(MMRSSE(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('Mean of Mean Relative SSE of three fields for different subspace dimension  -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('relative error')

%----------------------------
[num_temp,dim_temp]=size(MRMSSE_U);
figure
hold on
for i=1:num_temp
    semilogy(MRMSSE_U(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('MRMSSE of U field for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('relative error')

%----------------------------
[num_temp,dim_temp]=size(MRMSSE_V);
figure
hold on
for i=1:num_temp
    semilogy(MRMSSE_V(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('MRMSSE of V field for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('relative error')

%----------------------------
[num_temp,dim_temp]=size(MRMSSE_P);
figure
hold on
for i=1:num_temp
    semilogy(MRMSSE_P(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('MRMSSE of P field for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('relative error')

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


%Actual field plotting 
index=101;

Re=X_star(index,1);
u_lid=X_star(index,2);

U_in=U_star(index,:);
V_in=V_star(index,:);
P_in=P_star(index,:);
% P_in=P_in-min(P_in);    %Correction


figure
PlotUVPv(U_in,V_in,P_in,Re,u_lid,option_Data)
title(sprintf('%s GPE Predicted field -%d Training points -Re=%3.2f -ulid= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))

U_in=U_starorig(index,:);
V_in=V_starorig(index,:);
P_in=P_starorig(index,:);
P_in=P_in-min(P_in);    %Correction


figure
PlotUVPv(U_in,V_in,P_in,Re,u_lid,option_Data)
title(sprintf('Actual field  -%d Training points -Re=%3.2f -ulid= %3.3f ',num_train,X_star(index,1),X_star(index,2)))    



