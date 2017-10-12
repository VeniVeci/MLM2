% Test_PreImage_kD
% Test pure preimage on kpca and diffusion maps

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
num_train=50;                 
num_test=100;

%--------------------------------------------------------------------------
dim_new=20;           % temp value.


Koptions.ker='gaussian';
Koptions.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
Koptions.new_dim=dim_new;
Koptions.FullRec=0;

Doptions.metric ='euclidean';
Doptions.kernel ='gaussian'; 
% Doptions.kpara = 10000;             
Doptions.kAuto=1;

Doptions.dim_new = dim_new;              
Doptions.t = 1;                     
Doptions.FullRec = 0;      

% Diffusion PreImage options
preoptions.type='Dw';  %'LSE' OR 'Dw'
preoptions.para=2;
preoptions.neighbor=10;
% preoptions.neighbor=10;


SSE=[];

%%  test

    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);
    
    
%     [num,dim]=size(Y);
%     X_star=zeros(num,dim);  %Assign memory
%     
%     % kPCA
%     [Z_K,model_K] = Kpca(Y,KpcaOptions);
% 
%     %% Main      
%     for index=1:num %selecting every point in original space "X" to reconstruct  
%         Z_starN=Z_K(index,:);  
%         Do=zeros(num,1);      %Assign memory
%         for i =1:num          % Distance_original space
%             Do(i,1)  = Distance_OriginalSpace(Z_starN,i,model_K);    %Finding the distances original       
%         end 
%         X_star(index,:)= Dist2pos(model_K.X,Do,preoptions);
%         Anal.Distance_star(index,:)=Do';
%     end       
%     
%     
    % Auto kparameter for kPCA
%     Distance =pdist2(Y,Y,'euclidean');
%     Koptions.arg=sum(Distance(:))/(num_train^2);
    
    [Z_K,model_K] = Kpca2(Y,Koptions);
    [Z_D,model_D] = DiffusionMaps(Y,Doptions);
    
    
for i=1:dim_new
    
        Y_D_star = DiffusionMaps_PreImage(Z_D(:,1:i+1),model_D,preoptions);
        Y_K_star = Kpca_PreImage(Z_K(:,1:i),model_K,preoptions);
    
        SquErr_D=(Y-Y_D_star).^2;
        SquErr_K=(Y-Y_K_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);
        
        SSE_D(:,i)=(sum(SquErr_D'))';
        SSE_K(:,i)=(sum(SquErr_K'))';
        
end
    
 figure(1)
boxplot(SquErr_D')
title('Square error for all grid point -DiffusionMaps')

figure(2)
boxplot(SquErr_K')
title('Square error for all grid point -kPCA')

figure(3)
boxplot(SSE_D)
title('Square sum error for different subspace dimension -DiffusionMaps')
set(gca,'yscale','log');

figure(4)
boxplot(SSE_K)
title('Square sum error for different subspace dimension -kPCA')
set(gca,'yscale','log');

mean_SE_D=mean(SSE_D);
mean_SE_K=mean(SSE_K);

figure(5)
semilogy(mean_SE_D,'d');
hold on 
semilogy(mean_SE_K,'*');
hold off 
title('Mean of Square sum error for different subspace dimension -DiffusionMaps')
legend('DiffusionMaps','kPCA')
   
    
    
    