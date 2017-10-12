% Test_DiffusionMapVsPca


clear
%% Dataset Parameter
Name_field='exp1';
num_train=80;           %30;80 for exp1      %40 is bad  %Fine:50  %Best:80
num_test=400;


%--------------------------------------------------------------------------
new_dim=1;          % temp value.

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
DiffusionMap.Reoptions.para=5;
%  DiffusionMap.Reoptions.neighbor=6;

%--------------------------------------------------------------------------

%% Resolution Parameter
new_dim_start=1;    
new_dim_step =2;
new_dim_end  =9;

%% Plot Parameter 



%% Main program
% -----------------Defining the input ------------------------------------
options.order='normal';
[X,U,V,P,options] = Exp2UVPv(Name_field);
[X,X_star] = Data_Sepera(X,num_train,num_test,options);
[U,U_starorig] = Data_Sepera(U,num_train,num_test,options);
[V,V_starorig] = Data_Sepera(V,num_train,num_test,options);
[P,P_starorig] = Data_Sepera(P,num_train,num_test,options);


% [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
% [np_train,Dim_X]=size(X);
% [np_test, Dim_Y]=size(Y_starorig);

% -----------------Test several times-------------------------------------
for new_dim=new_dim_start:new_dim_step:new_dim_end
    
    % Diffusion maps GPR
    DiffusionMap.options.dim_new=new_dim;     
%   options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      
%   [Y_star_DiffusionMap,~,t_DiffusionMap]=GPR_DiffusionMap(X,Y,X_star,DiffusionMap.options,DiffusionMap.Reoptions);

    [ DiffusionMap.options.kpara ] = DiffusionMap_AutoK(U,DiffusionMap.options);
    [U_star_DiffusionMap,~,~]=GPR_DiffusionMap(X,U,X_star,DiffusionMap.options,DiffusionMap.Reoptions);
    
    [ DiffusionMap.options.kpara ] = DiffusionMap_AutoK(V,DiffusionMap.options);
    [V_star_DiffusionMap,~,~]=GPR_DiffusionMap(X,V,X_star,DiffusionMap.options,DiffusionMap.Reoptions);
    
    [ DiffusionMap.options.kpara ] = DiffusionMap_AutoK(P,DiffusionMap.options);
    [P_star_DiffusionMap,~,~]=GPR_DiffusionMap(X,P,X_star,DiffusionMap.options,DiffusionMap.Reoptions);
    
    % LPCA GPR
    [U_star_PCA,~,~]=GPR_SVD(X,U,X_star,new_dim);
    [V_star_PCA,~,~]=GPR_SVD(X,V,X_star,new_dim);
    [P_star_PCA,~,~]=GPR_SVD(X,P,X_star,new_dim);
    
    % Recording Result
    TempField_diff2=(U_starorig-U_star_DiffusionMap).^2;
    SSEU(:,new_dim)=sum(TempField_diff2,2);
    TempField_diff2=(V_starorig-V_star_DiffusionMap).^2;
    SSEV(:,new_dim)=sum(TempField_diff2,2);
    TempField_diff2=(P_starorig-P_star_DiffusionMap).^2;
    SSEP(:,new_dim)=sum(TempField_diff2,2);
    
    TempField_diff2=(U_starorig-U_star_PCA).^2;
    SSEUpca(:,new_dim)=sum(TempField_diff2,2);
    TempField_diff2=(V_starorig-V_star_PCA).^2;
    SSEVpca(:,new_dim)=sum(TempField_diff2,2);
    TempField_diff2=(P_starorig-P_star_PCA).^2;
    SSEPpca(:,new_dim)=sum(TempField_diff2,2);
    
    
      
%     % Recording Result
%     SquErr=(Y_starorig-Y_star_DiffusionMap).^2;
%     means=mean(Y_starorig,2);
%     SSErr=sum(SquErr,2);
%     RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
%     RateErr=abs(Y_starorig-Y_star_DiffusionMap)./Y_starorig;
%     RateE-rr=mean(RateErr,2);
% 
% 
%     RecSSErr_DiffusionMap(:,new_dim)   =SSErr;
%     RecRateSsErr_DiffusionMap(:,new_dim)=RatSsErr;
%     RecRateErr_DiffusionMap(:,new_dim)=RateErr;
%     RecTime_DiffusionMap(:,new_dim)=t_DiffusionMap;
%         

    
    
end

%% Plotting
plot_cursor=1;

   
   figure(plot_cursor)
   boxplot(SSEU(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
   title(sprintf('Diff-GPR SSE of each point; U field'));
   plot_cursor=plot_cursor+1;
   
   figure(plot_cursor)
   boxplot(SSEV(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
   title(sprintf('Diff-GPR SSE of each point; V field'));
   plot_cursor=plot_cursor+1;
   
   figure(plot_cursor)
   boxplot(SSEP(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
   title(sprintf('Diff-GPR SSE of each point; P field'));
   plot_cursor=plot_cursor+1;
   
   figure(plot_cursor)
   boxplot(SSEUpca(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
   title(sprintf('Pca-GPR SSE of each point; U field'));
   plot_cursor=plot_cursor+1;
   
   figure(plot_cursor)
   boxplot(SSEVpca(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
   title(sprintf('Pca-GPR SSE of each point; V field'));
   plot_cursor=plot_cursor+1;
   
   figure(plot_cursor)
   boxplot(SSEPpca(:,new_dim_start:new_dim_step:new_dim_end), {new_dim_start:new_dim_step:new_dim_end});
   title(sprintf('Pca-GPR SSE of each point; P field'));
   plot_cursor=plot_cursor+1;
   
   
   %%% X-Y Way %%% 
%    index=rand(5,1);
   index = randi([1 num_test],25,1);
   
   figure(plot_cursor)
   f=plot(U_star_DiffusionMap(index,:)',U_starorig(index,:)')
   hold on 
   plot([-10:10],[-10:10],'k*');
   plot([-10:10],[-10:10],'k');
   hold off
   axis equal
   axis([-2 10 -2 10])
%    x_size = 10;
%    y_size = 10;
%    set(f, 'PaperUnits', 'inches','PaperPosition',[0 0 x_size y_size]);
   title(sprintf('Duff-GPR U_star-U_original'))
   plot_cursor=plot_cursor+1;      
   
   figure(plot_cursor)
   plot(U_star_PCA(index,:)',U_starorig(index,:)')
   title(sprintf('Lpca-GPR U_star-U_original'))
   hold on 
   plot([-10:10],[-10:10],'k*');
   plot([-10:10],[-10:10],'k');
   hold off
   axis equal
   axis([-2 10 -2 10])
   plot_cursor=plot_cursor+1;   
   
   %%%%%%%%%%%%%%
      %%% X-Y Way %%% 
%    index=rand(5,1);
   index = randi([1 num_test],25,1);
   
   figure(plot_cursor)
   f=plot(V_star_DiffusionMap(index,:)',V_starorig(index,:)')
   hold on 
   plot([-10:10],[-10:10],'k*');
   plot([-10:10],[-10:10],'k');
   hold off
   axis equal
   axis([-6 3 -6 3])
%    x_size = 10;
%    y_size = 10;
%    set(f, 'PaperUnits', 'inches','PaperPosition',[0 0 x_size y_size]);
   title(sprintf('Duff-GPR U_star-U_original'))
   plot_cursor=plot_cursor+1;      
   
   figure(plot_cursor)
   plot(V_star_PCA(index,:)',V_starorig(index,:)')
   title(sprintf('Lpca-GPR U_star-U_original'))
   hold on 
   plot([-10:10],[-10:10],'k*');
   plot([-10:10],[-10:10],'k');
   hold off
   axis equal
   axis([-6 3 -6 3])
   plot_cursor=plot_cursor+1;   
   
   %%%%%%%%%%%%%%
   
   index = randi([1 num_test],5,1);
   figure(plot_cursor)
   plot(U_starorig(index,:)')
   hold on
   plot(U_star_DiffusionMap(index,:)','--') 
   plot(U_star_PCA(index,:)','*') 
   hold off
   title(sprintf('Coordinate to line; U field; Solid line for real case, dash for Duff-GPR prediction; * for Pca-GPR'))
   plot_cursor=plot_cursor+1;      
   
   index = randi([1 num_test],5,1);
   figure(plot_cursor)
   plot(V_starorig(index,:)')
   hold on
   plot(V_star_DiffusionMap(index,:)','--') 
   plot(V_star_PCA(index,:)','*') 
   hold off
   title(sprintf('Coordinate to line; V field; Solid line for real case, dash for Duff-GPR prediction; * for Pca-GPR'))
   plot_cursor=plot_cursor+1;      
   
   index = randi([1 num_test],5,1);
   figure(plot_cursor)   
   plot(P_starorig(index,:)')
   hold on
   plot(P_star_DiffusionMap(index,:)','--') 
   plot(P_star_PCA(index,:)','*') 
   hold off
   title(sprintf('Coordinate to line; P field; Solid line for real case, dash for Duff-GPR prediction; * for Pca-GPR'))
   plot_cursor=plot_cursor+1;      
   
   
   
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


