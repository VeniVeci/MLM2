% Test_Terminal_v1m2
% v1m1  version 1 modification 2. Add DiffusionMap GPR

clear
%% Dataset Parameter
filename='exp1';
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
DiffusionMap.options.kpara = 0.0000001;   % 100000 for CSTR; 10000 for originalCon; 10000 for originalSuperConductingWire100by100; 1 for original_TinMeltingFront100by100; 0.0000001 for 100100PloeFlow;     
DiffusionMap.options.dim_new = new_dim;              
DiffusionMap.options.t = 1;                     
DiffusionMap.options.FullRec = 0;      

% Diffusion PreImage options
DiffusionMap.Reoptions.type='Dw';
DiffusionMap.Reoptions.para=10;
% DiffusionMap.Reoptions.neighbor=5;

%--------------------------------------------------------------------------

%% Resolution Parameter
new_dim_start=1;    
new_dim_step =2;
new_dim_end  =5;

%% Plot Parameter 
boxplot_RateSsErr.active=1;
boxplot_RecRateErr.active=0;


%% Main program
% -----------------Defining the input ------------------------------------
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

% -----------------Test several times-------------------------------------
for new_dim=new_dim_start:new_dim_step:new_dim_end
    
    if lpca.active==1
        [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);
        
        % Recording Result
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
    
    if kpca.active==1;
        options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
        [Y_star_kpca,~,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
        
        % Recording Result
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
        [Y_star_isomap,~,t_isomap]=GPR_Isomap2(X,Y,X_star,isomap.options,isomap.Reoptions);
        
        % Recording Result
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
    
     if DiffusionMap.active==1;
        DiffusionMap.options.dim_new=new_dim;     
%         options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      
        [Y_star_DiffusionMap,~,t_DiffusionMap]=GPR_DiffusionMap(X,Y,X_star,DiffusionMap.options,DiffusionMap.Reoptions);
        
        % Recording Result
        SquErr=(Y_starorig-Y_star_DiffusionMap).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
        RateErr=abs(Y_starorig-Y_star_DiffusionMap)./Y_starorig;
        RateErr=mean(RateErr,2);
        
        RecSSErr_DiffusionMap(:,new_dim)   =SSErr;
        RecRateSsErr_DiffusionMap(:,new_dim)=RatSsErr;
        RecRateErr_DiffusionMap(:,new_dim)=RateErr;
        RecTime_DiffusionMap(:,new_dim)=t_DiffusionMap;
        
    end   
    
    
end

%% Plotting
plot_cursor=1;

if boxplot_RateSsErr.active==1
   if lpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));
       plot_cursor=plot_cursor+1;
   end
   
   if kpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));
       plot_cursor=plot_cursor+1;
   end
   
   if isomap.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_isomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Square Sum Error Rate of Each pixel of ISOMAP-GPR'));
       plot_cursor=plot_cursor+1;
   end
   
   if DiffusionMap.active==1;
       figure(plot_cursor)
       boxplot(RecRateSsErr_DiffusionMap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of DiffusionMap-GPR'));
       plot_cursor=plot_cursor+1;
   end   
   
   
end

if boxplot_RecRateErr.active==1
   if lpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_svd(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of LPCA-GPR'));
       plot_cursor=plot_cursor+1;
   end
   
   if kpca.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_kpca(:,new_dim_start:new_dim_step:new_dim_end),{new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of KPCA-GPR'));
       plot_cursor=plot_cursor+1;
   end
   
   if isomap.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_isomap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of ISOMAP-GPR'));
       plot_cursor=plot_cursor+1;
   end
   
   if DiffusionMap.active==1;
       figure(plot_cursor)
       boxplot(RecRateErr_DiffusionMap(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
       title(sprintf('Average Error Rate of Each pixel of DiffusionMap-GPR'));
       plot_cursor=plot_cursor+1;
   end   
   
   
   
   % A boxplot example
%     data = rand(20,24)
%     month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
%     simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];
%     boxplot(data,{month,simobs},'colors',repmat('rb',1,12),'factorgap',[5 2],'labelverbosity','minor');
    
%         %Plot
%         figure(new_dim)
%         title(sprintf('NonRescaled New dim= %g', new_dim))
%         plot(Y_star_kpca,'-ob')
%         hold on
%         plot(Y_star_svd,'-+r')
%         plot(Y_starorig,'-dk')
%         hold off
%         legend('Y^*-KPCA.','Y^*--LPCA','Y^*--real.','Location','northeast')


%         Mass=reshape(Y_star_svd,50,50);
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
end


