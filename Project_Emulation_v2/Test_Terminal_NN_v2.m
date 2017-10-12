% Test_Terminal_NN_v2
% Test Terminal Neural Network Version 2.0
%
% Pcakage Require:
% KPCA; PCA; Isomap MK-II;Neural Network;DiffusionMaps
% Example:
%
% About: 
% Modification:
% WeiX, May 5th 2014, First Edition

clear
%% Dataset Parameter
index_dataset=8;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR
num_train=40;                 
num_test=100;

%% Method Parameter
% for the moment kpca and DiffusionMap can not work at the same time due to
% some "same name" function crashion.
lpca.active=1;
kpca.active=0;
isomap.active=0;
DiffusionMap.active=0;

%--------------------------------------------------------------------------
% ANN Parameter
kfold=5;            % Value of cross validation. 1 means no cross validation happen.

%--------------------------------------------------------------------------
new_dim=1;          % temp value.

%--------------------------------------------------------------------------
%Lpca
lpca.options =0;    % No options for lpca

%--------------------------------------------------------------------------
%Kpca options
kpca.options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
%    options = struct('ker','rbf','arg',50,'new_dim',new_dim);                  % other options
%    options = struct('ker','poly','arg',[2,0],'new_dim',new_dim);              % other options
%    options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);    % other options
% isomap.options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      

%--------------------------------------------------------------------------
% Isomap options
isomap.options.dim_new=new_dim;                      % New dimension
isomap.options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
isomap.options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
isomap.options.metric='euclidean';             % Method of measurement. Metric

% Isomap PreImage options
isomap.Reoptions.ReCoverNeighborType='k';% Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
isomap.Reoptions.ReCoverNeighborPara=10;       % Parameter of neighbor of new point
isomap.Reoptions.Recoverd2pType='Dw';          % Type of distance to coordinate method. Distance weight/Least square estimate
isomap.Reoptions.Recoverd2pPara=1;             % Parameter of distance to coordinate recover method

%--------------------------------------------------------------------------
% DiffusionMap options
DiffusionMap.options.metric ='euclidean';
DiffusionMap.options.kernel ='gaussian'; 
DiffusionMap.options.kpara = 100000;   % 100000 for CSTR; 10000 for originalCon; 10000 for originalSuperConductingWire100by100; 1 for original_TinMeltingFront100by100; 0.0000001 for 100100PloeFlow;     
DiffusionMap.options.dim_new = new_dim;              
DiffusionMap.options.t = 1;                     
DiffusionMap.options.FullRec = 0;      

% Diffusion PreImage options
DiffusionMap.Reoptions.type='Dw';
DiffusionMap.Reoptions.para=10;
% DiffusionMap.Reoptions.neighbor=5;


%% Resolution Parameter
new_dim_start=1;    
new_dim_step =1;
new_dim_end  =3;

%% Plot Parameter 
boxplot_RateSsErr.active=1;
boxplot_RecRateErr.active=0;
lineplot.active=1;


%% Main program
% -----------------Defining the input ------------------------------------
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

% X=X';
% Y=Y';
% X_star=X_star';
% -----------------Test several times-------------------------------------
for new_dim=new_dim_start:new_dim_step:new_dim_end
    
    if lpca.active==1
        
        
        [Z_SVD,Key] = DataReduc_SVD(Y,new_dim);
        [Z_star_SVD]= ANN_prediction_crossvalidation(X',Z_SVD',X_star',kfold);

        Z_star_SVD=Z_star_SVD';
        Y_star_SVD=Z_star_SVD*Key;
        
        % Recording Result     
        SquErr=(Y_starorig-Y_star_SVD).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_SVD)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_svd(:,new_dim)   =SSErr;
        RecRateSsErr_svd(:,new_dim)=RatSsErr;
        RecRateErr_svd(:,new_dim)=RateErr;
        RecTime_svd(:,new_dim)=t_svd;     
        
    end
    
    if kpca.active==1;
%         options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
        kpca.options.new_dim=new_dim;  
        model2 = Kpca(Y,options);

        Z_star_kpca= ANN_prediction_crossvalidation(X,model2.Z,X_star,kfold);

        Y_star_kpca = Kpca_PreImage(Z_star_kpca,model2);
        
        [Y_star_kpca,t_kpca]=ANN_KPCA(X,Y,X_star,kpca.options,kfold); %5000000
        
        % Recording Result
        Y_star_kpca=Y_star_kpca';
        
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
        
    end
    
    if isomap.active==1;
        isomap.options.dim_new=new_dim;     
%         options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      
        [Y_star_isomap,t_isomap]=ANN_Isomap(X,Y,X_star,isomap.options,isomap.preoptions,kfold);
        
        % Recording Result
        Y_star_isomap=Y_star_isomap';
        
        SquErr=(Y_starorig-Y_star_isomap).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_isomap)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_isomap(:,new_dim)   =SSErr;
        RecRateSsErr_isomap(:,new_dim)=RatSsErr;
        RecRateErr_isomap(:,new_dim)=RateErr;
        RecTime_isomap(:,new_dim)=t_isomap;
        
    end
end

%% Plotting
plot_cursor=1;

if boxplot_RateSsErr.active==1
   if lpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Square Sum Error Rate of Each pixel of LPCA-ANN'));
       plot_cursor=plot_cursor+1;
   end
   
   if kpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Square Sum Error Rate of Each pixel of KPCA-ANN'));
       plot_cursor=plot_cursor+1;
   end
   
   if isomap.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_isomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Square Sum Error Rate of Each pixel of ISOMAP-ANN'));
       plot_cursor=plot_cursor+1;
   end
end

if boxplot_RecRateErr.active==1
   if lpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of LPCA-ANN'));
       plot_cursor=plot_cursor+1;
   end
   
   if kpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of KPCA-ANN'));
       plot_cursor=plot_cursor+1;
   end
   
   if isomap.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_isomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of ISOMAP-ANN'));
       plot_cursor=plot_cursor+1;
   end
end
   
if lineplot_instance.active==1
   figure(plot_cursor)
   y_plot=Y_starorig(lineplot_instance.index,:)';
   plot(y_plot);
   hold on
   
   if lpca.active==1;
       y_plot=Y_star_svd(lineplot_instance.index,:)';
       plot(y_plot,'--');
   end
   
   if kpca.active==1;
       y_plot=Y_star_kpca(lineplot_instance.index,:)';
       plot(y_plot,'*');
   end
   
   if isomap.active==1;
      y_plot=Y_star_isomap(lineplot_instance.index,:)';
       plot(y_plot,'+');
   end
   
   hold off
   title(sprintf('Acutal Instance Line Plot'));
   plot_cursor=plot_cursor+1;
end

   
   
   
   
