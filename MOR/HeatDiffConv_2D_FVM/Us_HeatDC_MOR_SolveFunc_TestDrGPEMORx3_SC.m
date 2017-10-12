%% Us_HeatDC_MOR_SolveFunc_TestDrGPEMORx3_SC
%
% Modifications:
% 5-4-2016, WeiX, first edition 

clear
tic
%%---------------Setting parameter---------------------------------------
Num_Train=80;   
% Index_test=303;
Index_test=54;

Num_Snapshot=10;

n_Ubases=10;     %Number of POD bases
dim_new=6;      % For LPCA more than 10 require. The new dimension of reduced emulatior

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


%% ----------------Load dataset-------------------------------------------
load('HeatDC_DBC_ExpData13.mat')
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best
    
% X_star=[10,-10,5];
Paras.t_n=Paras.t_end/Paras.dt;
Index_snapshot=1:(Paras.t_n/Num_Snapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

Y=Y_Rec(:,Index_snapshot,:);
Y=reshape(Y,[],500); %500 for this dataset
Y=Y';

X_star=X(Index_test,:);
Y_starorig=Y(Index_test,:);

X=X(1:Num_Train,:);
Y=Y(1:Num_Train,:);

% New desinated test point
% X_star(1,1)=0.55;
% X_star(1,2)=0.45;

uParas.ampx=1;
uParas.ampy=1;
uParas.x0=X_star(1,1);
uParas.y0=X_star(1,2);
%% ---------------HDM----------------------------------------------------
%-Disable As full data is provided.

tic;   
uParas.x0=X_star(1,1);
uParas.y0=X_star(1,2);

[Y_Full,~]=Us_HeatDC_SolveFunc(Paras,uParas,X_int);

Time_HDM=toc;

% Y_Full=Y_Rec(:,:,Index_test);
% 
% Time_Full=Time_Rec(Index_test);

%% ---------------MOR Bases by perfect snapshot -----------------------------------
toc0=toc;


Y_orig_snaps =reshape(Y_starorig',Paras.Num_node_y*Paras.Num_node_x,[]);
[U_OrigSnapS,S,~]=svd(Y_orig_snaps,'econ');  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
U_OrigSnapS=U_OrigSnapS(:,1:n_Ubases);
eigenvalues_OrigSnapS=diag(S);

[Y_U_OrigSnapS,Time_U_OrigSnapS]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U_OrigSnapS);

Time_ROM_U_OrigSnapS=toc-toc0;

%% ---------------MOR Golboal Bases -----------------------------------
Y_Global_snaps=reshape(Y_Rec,Paras.Num_node_y*Paras.Num_node_x,[]);
% n_Ubases=5;
% Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point
[U_GlobalSnapS,~,~]=svd(Y_Global_snaps,'econ');  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
U_GlobalSnapS=U_GlobalSnapS(:,1:n_Ubases);
% eigenvalues_GlobalSnapS=diag(S);


[Y_U_GlobalSnapS,Time_ROM_U_GlobalSnapS]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U_GlobalSnapS);



%  FiledPlot (Paras,T_orig,T_U_orig)


%% -------------MOR Bases by GPE Predicted snapshot-------------------------
% n_Ubases=5;   %number of POD basics

[Y_GPE,Time_Emu]=Func_DrGPE(X(1:Num_Train,:),Y(1:Num_Train,:),X_star,options,preoptions);     

Y_GPESnapS=reshape(Y_GPE,Paras.Num_node_y*Paras.Num_node_x,[]);
[U_GPESnapS,~,~]=svd(Y_GPESnapS);  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
U_GPESnapS=U_GPESnapS(:,1:n_Ubases);

[Y_U_GPESnapS,Time_ROM_U_GPESnapS]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U_GPESnapS);



%% -------------Figure local basis-------------------------
% Figure local basis
for i=1:Num_Train    
    Yi=Y(i,:);
    Yi_SnapS=reshape(Yi',Paras.Num_node_y*Paras.Num_node_x,[]);
    [U_YiSnapS,~,~]=svd(Yi_SnapS,'econ');  % U*S*V'=Rec_X
    U_YiSnapS=U_YiSnapS(:,1:n_Ubases); 
    U_Y(i,:)=U_YiSnapS(:)';
end

%% -------------MOR Bases by GPE Prediction-------------------------
[U_GPE,Time_Emu2]=Func_DrGPE(X(1:Num_Train,:),U_Y(1:Num_Train,:),X_star,options,preoptions);     

U_GPEi=reshape(U_GPE',Paras.Num_node_y*Paras.Num_node_x,[]);
[Y_U_GPE,Time_ROM_U_GPE]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U_GPEi);


%% ------------- Plot -------------------------
dx=Paras.Lx/Paras.Num_node_x;
dy=Paras.Ly/Paras.Num_node_y;
[Xcor,Ycor] = meshgrid(0+dx/2:dx:Paras.Lx-dx/2,0+dy/2:dy:Paras.Lx-dy/2);    

Index_frame=1;
for i=1:Paras.t_end/Paras.dt+1
    
%     Tx=reshape(P(1,:),Num_node_y,Num_node_x);       % Rearrange coordinate
%     Ty=reshape(P(2,:),Num_node_y,Num_node_x);

    T1 =reshape(Y_Full(:,i),Paras.Num_node_y,Paras.Num_node_x);
    T2=reshape(Y_U_OrigSnapS(:,i),Paras.Num_node_y,Paras.Num_node_x);
    T3=reshape(Y_U_GPESnapS(:,i),Paras.Num_node_y,Paras.Num_node_x);
    T4=reshape(Y_U_GlobalSnapS(:,i),Paras.Num_node_y,Paras.Num_node_x);
    T5=reshape(Y_U_GPE(:,i),Paras.Num_node_y,Paras.Num_node_x);
    
    figure(1)
    subplot(1,5,1),contour(Xcor,Ycor,T1),
    title('HDM'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,2),contour(Xcor,Ycor,T2),
    title('Local base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,3),contour(Xcor,Ycor,T3),
    title('GPE SnapS base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,4),contour(Xcor,Ycor,T4),
    title('Global base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,5),contour(Xcor,Ycor,T5),
    title('GPE base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    
    figure(2)
    subplot(1,5,1),pcolor(Xcor,Ycor,T1),
    title('HDM'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,2),pcolor(Xcor,Ycor,T2),
    title('Local base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,3),pcolor(Xcor,Ycor,T3),
    title('GPE SnapS base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,4),pcolor(Xcor,Ycor,T4),
    title('Global base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(1,5,5),pcolor(Xcor,Ycor,T5),
    title('GPE base'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal


     F(Index_frame)=getframe;
     Index_frame=Index_frame+1;

end

i=100;
% i=100;
T1 =reshape(Y_Full(:,i),Paras.Num_node_y,Paras.Num_node_x);
T2=reshape(Y_U_OrigSnapS(:,i),Paras.Num_node_y,Paras.Num_node_x);
T3=reshape(Y_U_GPESnapS(:,i),Paras.Num_node_y,Paras.Num_node_x);
T4=reshape(Y_U_GlobalSnapS(:,i),Paras.Num_node_y,Paras.Num_node_x);
T5=reshape(Y_U_GPE(:,i),Paras.Num_node_y,Paras.Num_node_x);

figure
pcolor(Xcor,Ycor,T1),shading interp;
% pcolor(Xcor-0.01,Ycor-0.01,T1),shading interp;

colormap(gray) %Only gray scale

title('(b)')


ax = gca; 
ax.FontSize=30;
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'x_2'},'FontSize',30,'FontWeight','bold');
xlabel({'x_1 '},'FontSize',30,'FontWeight','bold');

set(ax, 'xlim', [0 1])
set(ax, 'ylim', [0 1])

ax.XTick = [0:0.2:1];
ax.YTick = [0:0.2:1];

% caxis([0 50])      %color bar
% hcb=colorbar;
% set(hcb, 'ylim', [0 10])
% set(hcb,'YTick',[0:2:10])


%Adjust windows size 
h = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 700 700];   %Set position and size 
% set(h,'EdgeColor','none')


%% Relative error hear map plot
i=100;
% relativeError=abs(Y_U_GPESnapS(:,i)-Y_Full(:,i))./Y_Full(:,i);
absError=abs(Y_U_GPESnapS(:,i)-Y_Full(:,i));

T=reshape(absError,Paras.Num_node_y,Paras.Num_node_x);

figure
h=pcolor(Xcor,Ycor,T),shading interp;
colorbar 
% pcolor(Xcor-0.01,Ycor-0.01,T1),shading interp;

colormap(gray) %Only gray scale
colormap(hot)

title('(b)')

ax = gca; 
ax.FontSize=30;
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'x_2'},'FontSize',30,'FontWeight','bold');
xlabel({'x_1 '},'FontSize',30,'FontWeight','bold');

set(ax, 'xlim', [0 1])
set(ax, 'ylim', [0 1])

ax.XTick = [0:0.2:1];
ax.YTick = [0:0.2:1];

% caxis([0 50])      %color bar
% hcb=colorbar;
% set(hcb, 'ylim', [0 10])
% set(hcb,'YTick',[0:2:10])


%Adjust windows size 
h = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 700 700];   %Set position and size 
% set(h,'EdgeColor','none')

