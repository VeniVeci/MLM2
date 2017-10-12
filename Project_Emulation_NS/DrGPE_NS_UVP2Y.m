% DrGPE_NS_UVP2Y
% dimension reduction GPE for Navier
% stock equation Lid driven cavity for the special format of dataset( e.g. using U,V,P as Y
% which are multi-output.)_For all variable fields (U,V,P)

%% Initialize
clear
close all 


%% Dataset Parameter
Name_field='exp1';

[X_Data,U,V,P,option_Data] = Exp2UVPv(Name_field);

% Test_split=201;
% num_test=100;

%Test data
% X_star=X(Test_split:Test_split+num_test-1,:);
% U_starorig=U(Test_split:Test_split+num_test-1,:);
% V_starorig=V(Test_split:Test_split+num_test-1,:);
% P_starorig=P(Test_split:Test_split+num_test-1,:);

% num_train_max=num_train_ref(end);
%Train data

% U=U(1:num_train_max,:);
% V=V(1:num_train_max,:);
% P=P(1:num_train_max,:);

[~,Dim_U]=size(U);
[~,Dim_V]=size(V);
[~,Dim_P]=size(P);

Y_Data=[U,V,P];  

%% Data Generation
num_train=120;
num_test=300;
Test_split=201;

dim_new=5;   

% Test data
X_star=X_Data(Test_split:Test_split+num_test-1,:);
U_starorig=U(Test_split:Test_split+num_test-1,:);
V_starorig=V(Test_split:Test_split+num_test-1,:);
P_starorig=P(Test_split:Test_split+num_test-1,:);
Y_starorig=Y_Data(Test_split:Test_split+num_test-1,:);

% Train data
X=X_Data(1:num_train,:);
U=U(1:num_train,:);
V=V(1:num_train,:);
P=P(1:num_train,:);
Y=Y_Data(1:num_train,:);

%% DR method and parameters
DrMethod='kPCA';
% DrMethod='DiffusionMaps';

switch DrMethod
    
    case 'kPCA'
        options.ker='gaussian';
        options.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
        % Koptions.arg=500;    % 500 For index_dataset=8
        options.new_dim=dim_new;
        options.FullRec=0;       
    case 'DiffusionMaps'
        options.metric ='euclidean';
        options.kernel ='gaussian'; 
        % Doptions.kpara = 10000;             
        options.kAuto=1;
        options.dim_new = dim_new;              
        options.t = 1;                     
        options.FullRec = 0;      
    case 'Isomaps'
        
    otherwise 
        error('No such DR method')
end

% PreImage options----------------------
preoptions.type='Exp';  %'LSE' OR 'Dw'
% preoptions.para=2;
% preoptions.neighbor=10;



%% GPE structure and parameters

Dim_X=1; % only a temporary value as reminder. Would be detect later by the script.
covfunc = {@covSum,{@covSEard,@covNoise}}; 
hyp.cov = [zeros(Dim_X+1,1);0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];

InfMethod=@infExact;
% InfMethod=@infMCMC;


%% Main
% Initialize
  
% Dimension reduction

switch DrMethod
    case 'kPCA'
        % Auto select kernel parameter
        Distance =pdist2(Y,Y,'euclidean');
        options.arg=sum(Distance(:).^2)/(num_train^2);   
        options.arg=sqrt(options.arg/2);
        [Z,model] = Kpca(Y,options);    

    case 'DiffusionMaps'
        [Z,model] = DiffusionMaps(Y,options);

    case 'Isomaps'

    otherwise 
    error('No such DR method')
end
    
    
    [num_Z,dim_Z]=size(Z);    
%     Z_Rec(1:num_Z,1:dim_Z,i)=Z;

    % Assumened uncorrelated MISO GP on Z
    h2 = waitbar(0,'Loop2');
    for j=1:dim_Z
        
        %Initialize GPE 
        %Clean memory. !IMPORTANT! Change it with structure. 
        
%         [num_X,dim_X]=size(X);
%         hyp.cov = ones(dim_X+1,1);
%         sn = 0.1;
%         hyp.lik = log(sn);
%         hyp.mean=[];     
        
        [num_X,dim_X]=size(X);
        hyp.cov = [zeros(dim_X+1,1);0];
        hyp.lik = log(sn);
        hyp.mean=[];

        
        hyp = minimize(hyp, @gp, -100, InfMethod, meanfunc, covfunc, likfunc, X, Z(:,j));
%             hyp = minimize(hyp, @gp, -100, @infMCMC,  meanfunc, covfunc, likfunc, X, Z(:,i));
%             exp(hyp.lik)
%             nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i))
%             hypArray(i,j)=hyp;
        hyp_Rec(j)=hyp;
        [m(:,j) s(:,j)] = gp(hyp, InfMethod, meanfunc, covfunc, likfunc, X, Z(:,j), X_star);  
          
        waitbar(j/dim_Z,h2,'Loop2')
    
    end
    close(h2) 
    
    
    SSE=zeros(num_test,dim_new);
    
for k=dim_new:dim_new
        
        switch DrMethod    
            case 'kPCA'
                Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);                            
            case 'DiffusionMaps'
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
        
end
    
    figure
    boxplot(SMRSSE)
    title(sprintf('Boxplot of FSMRSSE -%s -%d Training points', DrMethod, num_train))
    set(gca,'yscale','log');
    xlabel('Dimensions')
    ylabel('SSE')
    
    
  
   %% 
   % Actual field plotting 
    index=108;
    
    Re=X_star(index,1);
    u_lid=X_star(index,2);
    
    U_in=U_star(index,:);
    V_in=V_star(index,:);
    P_in=P_star(index,:);
    
    figure
    PlotUVPv(U_in,V_in,P_in,Re,u_lid,option_Data)
    title(sprintf('%s GPE Predicted field -%d Training points -Re=%3.2f -ulid= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))
    
    U_in=U_starorig(index,:);
    V_in=V_starorig(index,:);
    P_in=P_starorig(index,:);
     
    figure
    PlotUVPv(U_in,V_in,P_in,Re,u_lid,option_Data)
    title(sprintf('Actual field  -%d Training points -Re=%3.2f -ulid= %3.3f ',num_train,X_star(index,1),X_star(index,2)))
    
    
    
    
    
    
%      VarRelaError=abs(s./m);
%      sumVarRelaError(i,:)=sum(VarRelaError);
%     figure(10+i)
%     boxplot(VarRelaError)
%     title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
%     set(gca,'yscale','log');
%     xlabel('Index of principal direction')
%     ylabel('Predicted variance ratio')
    



           
%     figure
%     boxplot(SSE)
%     title(sprintf('Boxplot of Square sum error for different subspace dimension -%s -%d Training points', DrMethod,num_train))
%     set(gca,'yscale','log');
%     xlabel('Dimension of manifold')
%     ylabel('Square sum error')

    
    
    
%     mean_SE(i,:)=mean(SSE);



% str=func2str(InfMethod);
% filename = sprintf('Data%d%s%s%s',index_dataset,DrMethod,str);
% save(filename)


