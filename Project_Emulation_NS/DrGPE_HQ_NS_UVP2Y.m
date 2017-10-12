% DrGPE_Dev_NS_UVP2Y
% A HQ terminal for dimension reduction GPE _Developer version _Special design for Navier
% stock equation for the special format of dataset( e.g. using U,V,P as Y
% which are multi-output.)_For all variable fields (U,V,P) combined
% toghther as a long vector 

%% Initialize
clear
close all 


%% Dataset Parameter


% [X,U,V,P,option_Data] = Exp2UVPv(Name_field);
% load('exp1.mat');

% clearvars -except option_Data
Name_field='exp1';
[X,U,V,P,option_Data] = Exp2UVPv(Name_field);

% option_Data.dt=dt;
% option_Data.lx=lx;
% option_Data.ly=ly;
% option_Data.nx=nx;
% option_Data.ny=ny;

%% Data Generation
num_train=120;
num_test=10;
Test_split=201;

dim_new=10;          

% Test data
X_star=X(Test_split:Test_split+num_test-1,:);
U_starorig=U(Test_split:Test_split+num_test-1,:);
V_starorig=V(Test_split:Test_split+num_test-1,:);
P_starorig=P(Test_split:Test_split+num_test-1,:);

% Train data
X=X(1:num_train,:);
U=U(1:num_train,:);
V=V(1:num_train,:);
P=P(1:num_train,:);

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
preoptions.neighbor=10;



%% GPE structure and parameters
% % Structure 1
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
% likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);


% %Structure 2
% %Isotropic covariance function
% % covfunc = @covSEiso; 
% % hyp.cov = [0; 0];
% 
% covfunc =@covSEard;
% hyp.cov = [];       % will be updated later 
% 
% likfunc = @likGauss; 
% sn = 0.1;
% hyp.lik = log(sn);
% 
% meanfunc=[];
% hyp.mean=[];


%Structure 3
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
[~,Dim_U]=size(U);
[~,Dim_V]=size(V);
[~,Dim_P]=size(P);

Y=[U,V,P];  

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
%Clean memory

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
        hyp = minimize(hyp, @gp, -100, InfMethod, meanfunc, covfunc, likfunc, X, Y(:,j));
        [m(:,j) s(:,j)] = gp(hyp, InfMethod, meanfunc, covfunc, likfunc, X, Y(:,j), X_star);          
        
%       waitbar(((i-1)*dim_new+j)/(length(num_train_ref)*dim_new),'Training System working very hard')
    
        waitbar(j/dim_Z,h2,'Loop2')
    
    end
    close(h2) 
    
    
%     SSE=zeros(num_test,dim_new);
    
    
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

        U_star=Y_star(:,1:Dim_U);
        V_star=Y_star(:,1+Dim_U:Dim_U+Dim_V);
        P_star=Y_star(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
        
        
        SquEr_U=(U_starorig-U_star).^2;
        SquEr_V=(V_starorig-V_star).^2;
        SquEr_P=(P_starorig-P_star).^2;

        SSE_U(:,k)=sum(SquEr_U');
        SSE_V(:,k)=sum(SquEr_V');
        SSE_P(:,k)=sum(SquEr_P');
        
    end
    
   %% 
   % Actual field plotting 
    index=4;
    
    Re=X(index,1);
    u_lid=X(index,2);
    
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


