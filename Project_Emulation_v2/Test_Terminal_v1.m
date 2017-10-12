% Test_Terminal_v1

clear
%% Dataset Parameter
index_dataset=5;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
num_train=100;                 
num_test=100;

%% Method Parameter
lpca.active=1;
kpca.active=0;
isomap.active=1;

% !! WARMING: change of the parameter(options) also needed to be made in "Main Prggram" if modification is required !!
new_dim=0;          % temp value.
lpca.options =0;    % No options for lpca
kpca.options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
%    options = struct('ker','rbf','arg',50,'new_dim',new_dim);                  % other options
%    options = struct('ker','poly','arg',[2,0],'new_dim',new_dim);              % other options
%    options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);    % other options
isomap.options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      

%% Resolution Parameter
new_dim_start=1;    
new_dim_step =1;
new_dim_end  =3;

%% Plot Parameter 
boxplot_RateSsErr.active=1;
boxplot_RecRateErr.active=1;


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
        [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
        
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
        options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      
        [Y_star_isomap,Yvar_star_isomap,t_isomap]=GPR_Isomap(X,Y,X_star,options);
        
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


