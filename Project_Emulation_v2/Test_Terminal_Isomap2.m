% Test_Terminal_Isomap2

% This terminal focus on comparation between 'k' Isomap and 'epsilon'
% Isomap.

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
            
num_train=40;     % Should be less than 200. num_train<=200            
num_test=100;

% num_train_start=30;    
% num_train_step =10;
% num_train_end  =100;


%% Method Parameter

%% Resolution Parameter
% new_dim_start=1;    
% new_dim_step =1;
% new_dim_end  =3;

%% Plot Parameter 


%% Main program

Dim_new_array=[3:2:3];

% -----------------Test several times-------------------------------------

% -----------------Defining the input ------------------------------------
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);
    

CircleIndex=1;
for new_dim = Dim_new_array  
    
    % ---------------- Isomap-k -----------------------------
    % parameters
%     Para_k_array=[10:10:50];
    Para_k_array=10;
%     new_dim=1;

    RecIndex=1;
    for para=Para_k_array

        isomap.options.dim_new=new_dim;                      % New dimension
        isomap.options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
        isomap.options.neighborPara=para;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
        isomap.options.metric='euclidean';             % Method of measurement. Metric

        % Isomap PreImage options
        isomap.Reoptions.ReCoverNeighborType='k';      % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
        isomap.Reoptions.ReCoverNeighborPara=para;       % Parameter of neighbor of new point
        isomap.Reoptions.Recoverd2pType='Dw';          % Type of distance to coordinate method. Distance weight/Least square estimate
        isomap.Reoptions.Recoverd2pPara=2;             % Parameter of distance to coordinate recover method


        [Y_star_isomapk,~,t_isomap]=GPR_Isomap2(X,Y,X_star,isomap.options,isomap.Reoptions);

        % Recording Result
        SquErr=(Y_starorig-Y_star_isomapk).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))

        RecSSErr_isomapk(:,RecIndex)   =SSErr;
        RecRateSsErr_isomapk(:,RecIndex)=RatSsErr;
    %     RecTime_isomapk(:,RecIndex)=t_isomap;      

        RecIndex=RecIndex+1;

    end
    
%     Index_LinePlot=randi([0,num_test],Num_LinePlot,1);
%     figure(2)
%     plot(Y_starorig(Index_LinePlot,:)')
%     hold on 
%     plot(Y_star_isomap(Index_LinePlot,:)','--')
%     hold off
%     title(sprintf('Line Plot of Graph of ISOMAP-k-GPR'));
    
    
    %--------------------------------------------------------

    % ---------------- Isomap-e -----------------------------
    % code helps finding range of distance
%     d=pdist2(Y,Y);
%     d=d(:);
%     figure(2);
%     hist(d(:,1));    
    % -----------
    
%     Para_e_array=[400:100:1000];
    Para_e_array=500;
%     new_dim=2;

    RecIndex=1;
    for para=Para_e_array

        isomap.options.dim_new=new_dim;                      % New dimension
        isomap.options.neighborType='epsilon';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
        isomap.options.neighborPara=para;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
        isomap.options.metric='euclidean';             % Method of measurement. Metric

        % Isomap PreImage options
        isomap.Reoptions.ReCoverNeighborType='epsilon';      % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
        isomap.Reoptions.ReCoverNeighborPara=para;       % Parameter of neighbor of new point
        isomap.Reoptions.Recoverd2pType='Dw';          % Type of distance to coordinate method. Distance weight/Least square estimate
        isomap.Reoptions.Recoverd2pPara=1;             % Parameter of distance to coordinate recover method


        [Y_star_isomape,~,t_isomap]=GPR_Isomap2(X,Y,X_star,isomap.options,isomap.Reoptions);

        % Recording Result
        SquErr=(Y_starorig-Y_star_isomape).^2;
        means=mean(Y_starorig,2);
        SSErr=sum(SquErr,2);
        RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))

        RecSSErr_isomape(:,RecIndex)   =SSErr;
        RecRateSsErr_isomape(:,RecIndex)=RatSsErr;
    %     RecTime_isomapk(:,RecIndex)=t_isomap;      

        RecIndex=RecIndex+1;
        
    end    
    
%     Index_LinePlot=randi([0,num_test],Num_LinePlot,1);
%     figure(3)
%     plot(Y_starorig(Index_LinePlot,:)')
%     hold on 
%     plot(Y_star_isomap(Index_LinePlot,:)','--')
%     hold off
%     title(sprintf('Line Plot of Graph of ISOMAP-e-GPR'));
    
    %--------------------------------------------------------    



    %% -------------------------Plotting -----------------------------
    figure(1)
    
    subplot(length(Dim_new_array),2,CircleIndex)
    boxplot(RecRateSsErr_isomapk, Para_k_array);
    title(sprintf('Square Sum Error Rate of Each pixel of ISOMAP-k-GPR'));
    ylabel('Square Sum Error Rate');xlabel('parameter');

    subplot(length(Dim_new_array),2,CircleIndex+1)
    boxplot(RecRateSsErr_isomape, Para_e_array);
    title(sprintf('Square Sum Error Rate of Each pixel of ISOMAP-e-GPR'));
    ylabel('Square Sum Error Rate');xlabel('parameter');

    
    Num_LinePlot=4;
    Index_LinePlot=randi([1,num_test],Num_LinePlot,1);
    figure(CircleIndex+1)
    plot(Y_starorig(Index_LinePlot,:)')
    hold on 
    plot(Y_star_isomapk(Index_LinePlot,:)','--')
    plot(Y_star_isomape(Index_LinePlot,:)','*')
    hold off
    title(sprintf('Line Plot of Graph of ISOMAP-GPR'));
    
    
    
    
    
    CircleIndex=CircleIndex+2;
    
   
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


%         Mass=reshape(Y_star_svd(1,:),100,100);
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



