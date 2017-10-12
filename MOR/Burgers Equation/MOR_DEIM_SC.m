%% MOR_DEIM_SC
%  version 2: using different dataset formate.
% Modifications:
% 7-2-2017, WeiX, first edition 

clear

%%---------------Setting parameter---------------------------------------
Num_Train=180;   
Num_Test=300;
Test_StartIndex=201;

Num_Snapshot=200;
Num_DEIMSnapshot=200;

n_Ubases=25;     %Number of POD bases
n_DEIM=30;

dim_new=10;      % For LPCA more than 10 require. The new dimension of reduced emulatior


%% ----------------Load dataset-------------------------------------------
% load('ExpDataV2_27.mat') 
load('ExpDataV2_24.mat') %CASE 2 g(x)=0.02*exp(x) para=[10,1000]x[2,5]x[0.2]
% load('ExpDataV2_9.mat') %Different data set has different initial and source term!!! Need change 'Burger1D_FEM_DBC_MOR_DEIM_SolverF' function  %CASE 1


%% DR method and parameters
     
DrMethod='kPCA';
% DrMethod='DiffusionMaps';
% DrMethod='Isomaps';
% DrMethod='PCA';

switch DrMethod
    
    case 'kPCA'
        options.ker='gaussian';   
        options.new_dim=dim_new;
        options.FullRec=0;       
        options.arg=1000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
        options.kAuto=1;
        
    case 'DiffusionMaps'
        options.metric ='euclidean';
        options.kernel ='gaussian'; 
        options.dim_new = dim_new;              
        options.t = 1;                     
        options.FullRec = 0;      
        % Doptions.kpara = 10000;             
        options.kAuto=1;
        options.Ztype=0;    %Output type. With/without 1st component
        
    case 'Isomaps'
        options.dim_new=dim_new;                % New dimension
        options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
        options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
        options.metric='euclidean';             % Method of measurement. Metric

        %Isomap PreImage options
        preoptions.ReCoverNeighborType='k';     % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
        preoptions.ReCoverNeighborPara=10;      % Parameter of neighbor of new point
        
    case 'PCA'
        options=[];
        
    otherwise 
        error('No such DR method')
end


options.DrMethod=DrMethod;
% PreImage options----------------------
preoptions.type='Exp';  %'LSE', 'Dw' OR 'Exp'
preoptions.neighbor=10;
% preoptions.type='LpcaI';
% preoptions.dim_new=10; % Use to stable the result but sacrefy accuracy
% preoptions.InMethod='ANN';

%--------------
Index_snapshot=1:(Paras.t_n/Num_Snapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

Y=Y_Rec(:,Index_snapshot,:);
Y=reshape(Y,[],500); %500 for this dataset
Y=Y';
Time_HDM=Time;

X_star=X(Test_StartIndex:Test_StartIndex+Num_Test-1,:);
Y_starorig=Y(Test_StartIndex:Test_StartIndex+Num_Test-1,:);

X=X(1:Num_Train,:);
Y=Y(1:Num_Train,:);


%---------------
Index_DEIMsnapshot=1:(Paras.t_n/Num_DEIMSnapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

Y_DEIM=Y_Rec(:,Index_DEIMsnapshot,:);
Y_DEIM=reshape(Y_DEIM,[],500); %500 for this dataset
Y_DEIM=Y_DEIM';

Y_DEIM=Y_DEIM.^2;

Y_DEIM_starorig=Y_DEIM(Test_StartIndex:Test_StartIndex+Num_Test-1,:);


Y_DEIM=Y_DEIM(1:Num_Train,:);



%% -------------MOR Bases by GPE Predicted snapshot-------------------------
% n_Ubases=5;   %number of POD basics
[Y_GPE,Time_Emu]=Func_DrGPE(X(1:Num_Train,:),Y(1:Num_Train,:),X_star,options,preoptions);     

[Y_DEIM_GPE,Time_Emu2]=Func_DrGPE(X(1:Num_Train,:),Y_DEIM(1:Num_Train,:),X_star,options,preoptions);     

%     %Global DEIM
%     Ysqr_GPESnapS=reshape(Y_DEIM',Paras.n+1,[]);
% 
%     [U_DEIM,S_DEIM,~]=svd(Ysqr_GPESnapS(2:end-1,:),'econ'); 
% 
%     [~,U_DEIM,P] = DEIM(U_DEIM);



h = waitbar(0,'Test MOR Bases by GPE snapshot');
for i =1:Num_Test
        
            
        Paras.Re=X_star(i,1);
        Paras.u0a=X_star(i,2);
        Paras.u0b=X_star(i,3);
       
        Y_GPESnapS=reshape(Y_GPE(i,:)',Paras.n+1,[]);
        Ysqr_GPESnapS=reshape(Y_DEIM_GPE(i,:)',Paras.n+1,[]);
        
        [U_GPESnapS,~,~]=svd(Y_GPESnapS(2:end-1,:));  % U*S*V'=Rec_X   
        
        [U_DEIM,S_DEIM,~]=svd(Ysqr_GPESnapS(2:end-1,:)); 
        
        [~,U_DEIM,P] = DEIM(U_DEIM);

       
        
        for j=[25,35]
        
            U_GPESnapSj=U_GPESnapS(:,1:j);
%             U_DEIMj=U_DEIM(:,1:j);

%             n_DEIM=j; % n_DEIM =n_POD
                        
            U_DEIMj=U_DEIM(:,1:n_DEIM);                      
            Pj=P(:,1:n_DEIM);            %same number as n_Ubases
                 
            [Y_U_GPESnapS,~,~]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U_GPESnapSj,U_DEIMj,Pj);
        
            Y_U_GPESnapS_Rec(:,:,i)=Y_U_GPESnapS;
            Y_U_GPESnapS_Rec2(:,:,i,j)=Y_U_GPESnapS;
            
        
            SSE_dx_U_GPESnapS(i,:)=sum((Y_U_GPESnapS_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1);   %Square sum error; integral on dx
            SSE_dxdt_U_GPESnapS(i,j)=sum(SSE_dx_U_GPESnapS(i,:),2);                              %Square sum error; integral on dx & dt
            RE_U_GPESnapS(i,j)=mean(sqrt(SSE_dx_U_GPESnapS(i,:) ./ sum(Y_Rec(:,:,Test_StartIndex+i-1).^2,1)));
            
%             RE_dxdt_U_GPESnapS(i,j)
            
        
        end
    waitbar(i/Num_Test);
   
end
close(h);


%Take 5 largest error
RE_U_GPESnapS_Sort=sort(RE_U_GPESnapS);
RE_U_GPESnapS_Sort=RE_U_GPESnapS_Sort(1:end-5,:);



figure
boxplot(RE_U_GPESnapS_Sort)
title(sprintf('L^2 Error Boxplot. %s+DEIM, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i,Num_{SnapshotDEIM}=%0i Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,Num_DEIMSnapshot,dim_new))
set(gca,'yscale','log');
% ylabel({'Relative error'},'FontSize',30,'FontName','Times New Roman');
% xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman');


figure
boxplot(RE_U_GPESnapS(:,1:1:n_Ubases))
xlabel('POD basis dimension');
ylabel('Relative error');
title('(a)');

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
%     uB=0;
%     dB=-5;
% 
%     ax.YLim=[10^dB,10^uB];
%     ytick=logspace(dB,uB,uB-dB+1);

%     ytick=[,logspace(dB+1,uB,uB-dB+1-1)];

%     ax.YTick=ytick(1:1:end);
%     % ax.YTickLabel=ytick(1:2:end);
% 
%     % ax.XTick=5:10:n_Ubases-5;
%     ax.XTick=5:5:30;
%     ax.XTickLabel=5:5:n_Ubases;
% 
%     ax.YMinorTick = 'on';
%     ax.YGrid = 'on';




% xtickLabel=5:10:n_Ubases-5;
% xtickLabel=blanks(12);

% xtickLabel = cell(12,1);
% xtickLabel(:) = {''};
% xtickLabel(1:2:11)=5:10:n_Ubases-5;
% ax.XTickLabel=xtickLabel;




% ax.XTickLabel={"xtickLabel"};

% ax.XTickLabel=5:10:n_Ubases-5;


% ax.XLabel.L=5:5:n_Ubases;
% xt=get(gca,'xtick'); % get the current ones
% xt(1:2:end)=[]; % wipe out every other one
% set(gca,'xtick',xt) % and replace w/ updated

%% Analysis for representative case 
mColumn=25;
data=RE_U_GPESnapS(:,mColumn);

meanData=mean(data);
minData=min(data);
maxData=max(data);
quantileData=quantile(data,[.25 .5 .75]);


CompareValue=1.5*quantile(data,[.75]); %Upper whisker
% CompareValue=maxData;
CompareValue=meanData;

sortData=(data-CompareValue).^2;
[sortData,sortData_Index]=sort(sortData);
SortX_star=X_star(sortData_Index,:);

% Case plot

Test_index=sortData_Index(1);
y_GPE_DEIM=Y_U_GPESnapS_Rec2(:,:,Test_index,mColumn);

y_orig=Y_Rec(:,:,Test_StartIndex+Test_index-1);

% y_GPE_DEIM=Y_U_GPESnapS_Rec(:,:,Test_index);

h=1/Paras.n;      % space step size
x = 0:h:1;  % coordinate sequence

number_show=5;
figure 
hold on 

index=1:Paras.t_n/number_show:Paras.t_n+1;
index=[0,10,20,30,40,50,100,150,200]+1;

for i = index
    plot(x,y_orig(:,i),'k--','LineWidth',1.5)
    plot(x,y_GPE_DEIM(:,i),'k-','LineWidth',1.5)

end 
box on
hold off

ax = gca; 
ax.YLabel.Position=[-0.05,0.05,-1];
ax.YLabel.Position=[-0.07,0.1,-1];



ax = gca; 
ax.YTick=[-0.3:0.3:0.6];
ax.XTick=[0:0.2:1];
ylim([-0.3,0.6]);




ax = gca; 
ax.FontSize=30;
ax.FontName='Times New Roman';
ax.FontWeight='normal';
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'u(t, x)'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
xlabel({'x'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
title('(r=25 mean)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')
% 
% title(sprintf('Actual plot case %0i for %s Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{DimNew}=%i POD=%i', Test_index, DrMethod,Num_Train,Num_Test,Num_Snapshot,Num_DEIMSnapshot,dim_new,n_Ubases))



%     ax = gca; 
%     ax.YLabel.Position=[-0.05,0.05,-1];
%     ax.YLabel.Position=[-0.07,0.1,-1];
% 
% 
% 
%     ax = gca; 
%     ax.YTick=[-0.3:0.3:0.6];
%     ax.XTick=[0:0.2:1];
%     ylim([-0.3,0.6]);
% 
%     ax = gca; 
%     ax.YTick=[0:0.1:0.3];
%     ax.XTick=[0:0.2:1];
%     ylim([0,0.3]);






