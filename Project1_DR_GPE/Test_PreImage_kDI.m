% Test_PreImage_kDI
% Test pure preimage on kpca and diffusion maps and Isomaps

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
num_train=100;                 
num_test=1;

%--------------------------------------------------------------------------
dim_new=20;           % temp value.

% kernel PCA options----------------------
Koptions.ker='gaussian';
Koptions.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
% Koptions.arg=500;    % 500 For index_dataset=8
Koptions.new_dim=dim_new;
Koptions.FullRec=0;


% Diffusion Maps options--------------------------------------
Doptions.metric ='euclidean';
Doptions.kernel ='gaussian'; 
% Doptions.kpara = 10000;             
Doptions.kAuto=1;
Doptions.dim_new = dim_new;              
Doptions.t = 1;                     
Doptions.FullRec = 0;      

% Isomaps options------------------------------------------------------
Ioptions.dim_new=dim_new;                      % New dimension
Ioptions.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
Ioptions.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
Ioptions.metric='euclidean';             % Method of measurement. Metric

% PreImage options------------------------------------------------------
preoptions.type='Dw';  %'LSE' OR 'Dw'
preoptions.para=2;
preoptions.neighbor=10;
% preoptions.neighbor=10;

%%  test

    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);
    
    
    [Z_K,model_K] = Kpca2(Y,Koptions);
    [Z_D,model_D] = DiffusionMaps(Y,Doptions);
    [Z_I,model_I] = Isomaps(Y,Ioptions);
    
    
    for i=1:dim_new
    
        % curropt data
%         Z_K=Z_K*1.01;
%         Z_D=Z_D*1.01;
%         Z_I=Z_I*1.01;
%         Y=Y*1.01;
%         Y_D_star = DiffusionMaps_PreImage(Z_D(:,1:i+1),model_D,preoptions);
%         Y_K_star = Kpca_PreImage(Z_K(:,1:i),model_K,preoptions);
%         Y_I_star = Isomaps_PreImage(Z_I(:,1:i),model_I,preoptions);
               
        for index =1:num_train
            Z_star=Z_I(index,1:i);
            model_temp=model_I;       % model_temp is model without the index point to avoide 100% accurate recovery.
            model_temp.Z=removerows(Z_I,'ind',index);  
            model_temp.X=removerows(Y,'ind',index);  
            Y_I_star(index,:) = Isomaps_PreImage(Z_star,model_temp,preoptions);
        end
        
       
        SquErr_D=(Y-Y_D_star).^2;
        SquErr_K=(Y-Y_K_star).^2;
        SquErr_I=(Y-Y_I_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);
        
        SSE_D(:,i)=(sum(SquErr_D'))';
        SSE_K(:,i)=(sum(SquErr_K'))';
        SSE_I(:,i)=(sum(SquErr_I'))';
        
end
    
SquErr_K=real(SquErr_K);
SSE_K=real(SSE_K);


figure(1)
boxplot(SquErr_D')
title('Square error for all grid point for each test point -DiffusionMaps')
set(gca,'yscale','log');

figure(2)
boxplot(SquErr_K')
title('Square error for all grid point for each test point -kPCA')
set(gca,'yscale','log');

figure(3)
boxplot(SquErr_I')
title('Square error for all grid point for each test point -Isomap')
set(gca,'yscale','log');



figure(4)
boxplot(SSE_D)
title('Square sum error for different subspace dimension -DiffusionMaps')
set(gca,'yscale','log');

figure(5)
boxplot(SSE_K)
title('Square sum error for different subspace dimension -kPCA')
set(gca,'yscale','log');

figure(6)
boxplot(SSE_I)
title('Square sum error for different subspace dimension -Isomaps')
set(gca,'yscale','log');


mean_SE_D=mean(SSE_D);
mean_SE_K=mean(SSE_K);
mean_SE_I=mean(SSE_I);


figure(7)
semilogy(mean_SE_D,'d');
hold on 
semilogy(mean_SE_K),'*';
semilogy(mean_SE_I,'o');
hold off 
title('Mean of Square sum error for different subspace dimension -DiffusionMaps')
legend('DiffusionMaps','kPCA','Isomaps')    
    
    