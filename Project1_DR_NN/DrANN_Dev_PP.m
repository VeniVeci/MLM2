% DrANN_Dev
% A HQ terminal for dimension reduction ANN with Pre-processing

%% Initialize
clear
close all 


%% Dataset Parameter
index_dataset=5;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR
% num_train=200;                 

% num_train_ref=[50:50:150]'; % must be in INCREASING order

num_train_ref=[40:40:200]'; % must be in INCREASING order
num_test=300;

dim_new=20;           

kfold=1;

%% DR method and parameters
DrMethod='kPCA';
DrMethod='DiffusionMaps';

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



%% NN structure and parameters


%% Main
% Initialize
Z_Rec=zeros(num_train_ref(end),dim_new+1,length(num_train_ref));    % +1 just for diffusion maps
h1 = waitbar(0,'Loop1');

for i=1:length(num_train_ref)
    
    waitbar(i/length(num_train_ref),h1)
    
    % Generate dataset
    num_train=num_train_ref(i);
    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

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
    Z_Rec(1:num_Z,1:dim_Z,i)=Z;
%     model_Rec(i)=model;
    
    % Assumened uncorrelated MISO GP on Z
    h2 = waitbar(0,'Loop2');
    
    

    
    for k=1:dim_new
          
        switch DrMethod    
            case 'kPCA'      
                
                [Z_train,key] = DataPP(Z(:,1:k));                  
                Z_star = ANN_prediction_crossvalidation(X',Z_train',X_star',kfold);
                Z_star=Z_star';
                Z_star = DataPP_Rocv(Z_star,key);
                Y_star = Kpca_PreImage(Z_star,model,preoptions);
                                
            case 'DiffusionMaps'             
                [Z_train,key] = DataPP(Z(:,2:k+1));               
                Z_star= ANN_prediction_crossvalidation(X',Z_train',X_star',kfold);
                Z_star=Z_star'; 
                Z_star = DataPP_Rocv(Z_star,key);
                Z_star=[repmat(Z(1,1),[num_test,1]),Z_star];           
                Y_star = DiffusionMaps_PreImage(Z_star,model,preoptions);
                
%             case 'DiffusionMaps'             
%                 [Z_train,key] = DataPP(Z(:,1:k+1));               
%                 Z_star= ANN_prediction_crossvalidation(X',Z_train',X_star',kfold);
%                 Z_star=Z_star'; 
%                 Z_star = DataPP_Rocv(Z_star,key);     
%                 Y_star = DiffusionMaps_PreImage(Z_star,model,preoptions);    
                      
            case 'Isomaps'

            otherwise 
                error('No such DR method')
        end

        SquErr=(Y_starorig-Y_star).^2;        
        MSS_orig=mean(Y_starorig.^2,2);
        
        SSE(:,k)=sum(SquErr,2);
        RSSE(:,k)=mean(SquErr,2)./MSS_orig;

        waitbar(k/dim_Z,h2,'Loop2')
    
    end
    close(h2) 
   
    
    figure(i)
    boxplot(RSSE)
    title(sprintf('Boxplot of RSSE for different subspace dimension -%s -%d Training points', DrMethod, num_train_ref(i)))
    set(gca,'yscale','log');
    xlabel('Dimensions')
    ylabel('RSSE')

     MRSSE(i,:)=mean(RSSE);
    
end

close(h1) 


markers = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
[num_temp,dim_temp]=size(MRSSE);
figure
hold on
for i=1:num_temp
    semilogy(MRSSE(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('Relative error (MRSSE=Mean SSE/Morig^2) for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('Relative error')

% str=func2str(InfMethod);
% filename = sprintf('Data%d%s%s%s',index_dataset,DrMethod,str);
% save(filename)


