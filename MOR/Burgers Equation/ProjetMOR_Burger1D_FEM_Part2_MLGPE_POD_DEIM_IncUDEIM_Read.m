%% ProjetMOR_Burger1D_FEM_Part2_MLGPE_POD_DEIM_IncUDEIM_Read
%  version 2: using different dataset formate.
% Modifications:
% 7-Oct-2016, WeiX, first edition 

clear



%% ----------------Load dataset-------------------------------------------
load('Bur_MOR_kGPE_DEIM_Train180Test300SS200DEIMSS200U30UDEIM5to60v73.mat') %7 is ok for HH, HH fails in 4.
% load('Bur_MOR_kGPE_DEIM_Train180Test300SS200DEIMSS200U50UDEIM5to60v73.mat') 


h = waitbar(0,'Test MOR Bases by GPE snapshot');
for i =1:Num_Test
        
  
        for j=nDEIM
            
            Y_U_GPESnapS=Y_U_GPESnapS_Rec2(:,:,i,j);
        
            SSE_dx_U_GPESnapS(i,:)=sum((Y_U_GPESnapS-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1);   %Square sum error; integral on dx
            
            SSE_dxdt_U_GPESnapS(i,j)=sum(SSE_dx_U_GPESnapS(i,:),2);                              %Square sum error; integral on dx & dt
            
%             RE_U_GPESnapS(i,j)=mean (sqrt(sum((Y_U_GPESnapS_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1) ./ sum(Y_Rec(:,:,Test_StartIndex+i-1).^2,1)) ); 
            
            RE_U_GPESnapS(i,j)=mean(sqrt(SSE_dx_U_GPESnapS(i,:) ./ sum(Y_Rec(:,:,Test_StartIndex+i-1).^2,1)));
        
        end
    waitbar(i/Num_Test);
   
end
close(h);

figure
boxplot(SSE_dxdt_U_GPESnapS)
title(sprintf('L^2 Error Boxplot. %s+DEIM, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i,Num_{SnapshotDEIM}=%0i Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,Num_DEIMSnapshot,dim_new))
set(gca,'yscale','log');
% ylabel({'Relative error'},'FontSize',30,'FontName','Times New Roman');
% xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman');

nDEIM=1:n_DEIM;

figure
boxplot(SSE_dxdt_U_GPESnapS(:,nDEIM))
xlabel('DEIM basis dimension');
ylabel('Square error');
title('(b)');

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 


figure
boxplot(RE_U_GPESnapS(:,nDEIM))
xlabel('DEIM basis dimension');
ylabel('Relative error');
title('(b)');

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=1;
dB=-4;

ax.YLim=[10^dB,10^uB];
ytick=logspace(dB,uB,uB-dB+1);

ax.YTick=ytick(1:2:end);
% ax.YTickLabel=ytick(1:2:end);

% ax.XTick=5:10:n_Ubases-5;
ax.XTick=1:2:13;
ax.XTickLabel=nDEIM;
% ax.XTickLabel=5:10:55;
ax = gca; 
ax.XTickLabel=5:10:55;

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
mColumn=15;
data=RE_U_GPESnapS(:,mColumn);

meanData=mean(data);
minData=min(data);
maxData=max(data);
quantileData=quantile(data,[.25 .5 .75]);


CompareValue=1.5*quantile(data,[.75]); %Upper whisker
%  CompareValue=meanData;
% CompareValue=maxData;


sortData=(data-CompareValue).^2;
[sortData,sortData_Index]=sort(sortData);
SortX_star=X_star(sortData_Index,:);

RE=data(sortData_Index(1))
xStar=X_star(sortData_Index(1),:)


%% Case plot

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


markerIndex=1:3:length(x);
for i = index
%     plot(x,y_orig(:,i),'k--','LineWidth',1.5)
%     plot(x,y_orig(:,i),'ko','MarkerFaceColor','k','LineWidth',1)
    plot(x(markerIndex),y_orig(markerIndex,i),'ko','MarkerSize',5)
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
title('(b)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')
% 
% title(sprintf('Actual plot case %0i for %s Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{DimNew}=%i POD=%i', Test_index, DrMethod,Num_Train,Num_Test,Num_Snapshot,Num_DEIMSnapshot,dim_new,n_Ubases))



ax = gca; 
ax.YLabel.Position=[-0.05,0.05,-1];
ax.YLabel.Position=[-0.07,0.1,-1];



ax = gca; 
ax.YTick=[-0.3:0.3:0.6];
ax.XTick=[0:0.2:1];
ylim([-0.3,0.6]);


ax = gca; 
ax.YTick=[0:0.1:0.2];
ax.XTick=[0:0.2:1];
ylim([0,0.2]);

%%

absError=abs(y_GPE_DEIM-y_orig);
% y_GPE_DEIM=Y_U_GPESnapS_Rec(:,:,Test_index);

h=1/Paras.n;      % space step size
x = 0:h:1;  % coordinate sequence

number_show=5;
figure 
hold on 

index=1:Paras.t_n/number_show:Paras.t_n+1;
index=[0,10,20,30,40,50,100,150,200]+1;


markerIndex=1:3:length(x);
for i = index
%     plot(x,y_orig(:,i),'k--','LineWidth',1.5)
%     plot(x,y_orig(:,i),'ko','MarkerFaceColor','k','LineWidth',1)

    plot(x,absError(:,i),'-','LineWidth',1.5)

end 
box on
hold off


%% show particular case
% Test_index=40;
% y_orig=Y_Rec(:,:,Test_StartIndex+Test_index-1);
% y_GPE_DEIM=Y_U_GPESnapS_Rec(:,:,Test_index);
% 
% h=1/Paras.n;      % space step size
% x = 0:h:1;  % coordinate sequence
% figure
% for i = 1:Paras.t_n+1
%     
%     subplot(2,1,1)
%     plot(x,y_orig(:,i))
% %     axis([0,1,0,1]);
%     title(sprintf('FEM Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))
%     
%     subplot(2,1,2)
%     plot(x,y_GPE_DEIM(:,i))
% %     axis([0,1,0,1]);
%     title(sprintf('MOR FEM Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))   
% 
%     F(i) = getframe;
%     
% end


%%
% Make sure the index end means the same frame!!!

% figure(1)
% L1=plot(Y_orig(:),'k-');
% hold on
% L2=plot(Y_star_snaps(:),'b--');
% L3=plot(Y_U_star(:),'r-.');
% legend([L1,L2,L3],'Original Y filed','GPR Predicted Y filed','GPR-MOR Y filed')
% hold off




%% Error Analysis
% tic;
% h = waitbar(0,'Error analysis');
% for i =1:Num_Test
%     
%     SSE_dx_U_OrigSnapS(i,:)=sum((Y_U_OrigSnapS_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1);   %Square sum error; integral on dx
%     SSE_dx_U_GlobalSnapS(i,:)=sum((Y_U_GlobalSnapS_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1);
%     SSE_dx_U_GPESnapS(i,:) =sum((Y_U_GPESnapS_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1); 
%     SSE_dx_Y_U_GPE(i,:)=sum((Y_U_GPE_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1);
%     
%     SSE_dxdt_U_OrigSnapS(i,:)=sum(SSE_dx_U_OrigSnapS(i,:),2);                              %Square sum error; integral on dx & dt
%     SSE_dxdt_U_GlobalSnapS(i,:)=sum(SSE_dx_U_GlobalSnapS(i,:),2);                              %Square sum error; integral on dx & dt
%     SSE_dxdt_U_GPESnapS(i,:)=sum(SSE_dx_U_GPESnapS(i,:),2);                              %Square sum error; integral on dx & dt
%     SSE_dxdt_Y_U_GPE(i,:)=sum(SSE_dx_Y_U_GPE(i,:),2);                              %Square sum error; integral on dx & dt
%     
%     waitbar(i/Num_Test);
% end
% close(h);
% toc

% SE_U_OrigSnapS=(Y_U_OrigSnapS-Y_Full).^2;
% SE_U_GPESnapS=(Y_U_GPESnapS-Y_Full).^2;
% SE_U_GlobalSnapS=(Y_U_GlobalSnapS-Y_Full).^2;
% SE_U_GPE=(Y_U_GPE-Y_Full).^2;
% 
% SE_Y_GPESnapS=(Y_GPESnapS-Y_orig_snaps).^2;
% SE_Y_GPESnapS =repmat(SE_Y_GPESnapS,Paras.t_n/Paras.n_snap,1);
% SE_Y_GPESnapS=reshape(SE_Y_GPESnapS(:),Paras.n+1,[]);

%% ---- plot------------------------------------------------------
% ----Box plot------------------------------------------------------------
% Err_Y_star_snaps=(Y_GPESnapS-Y_orig_snaps).^2;
% Err_Y_U_orig=(Y_U_OrigSnapS-Y_Full).^2;
% Err_Y_U_star=(Y_U_GPESnapS-Y_Full).^2;

% figure
% boxplot([SSE_dxdt_U_OrigSnapS,SSE_dxdt_U_GlobalSnapS,SSE_dxdt_U_GPESnapS,SSE_dxdt_Y_U_GPE],'labels',{'ROM ferfect base','ROM global base','ROM GPE Snapshot base','ROM GPE base'})
% title(sprintf('L^2 Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,n_Ubases,dim_new))
% hold off
% set(gca,'yscale','log');
% figure(3)
% boxplot(Err_Y_star_snaps(:),'labels',{'Emulation Solution'}) 
% title(sprintf('L^2 Error on each point at each time'))


% % ----L2 Error accumulation plot------------------------------------------------------------
% i = 1:Paras.t_n+1;
% x=(i-1)*(Paras.t_end/Paras.t_n);
% figure(4)
% 
% L1=plot(x,cumsum(sum(SE_U_OrigSnapS)),'b--');
% hold on
% L2=plot(x,cumsum(sum(SE_U_GPESnapS)),'r-.');
% L3=plot(x,cumsum(sum(SE_U_GlobalSnapS)),'g--');
% L4=plot(x,cumsum(sum(SE_U_GPE)),'y--');
% L5=plot(x,cumsum(sum(SE_Y_GPESnapS(:,1:Paras.t_n+1))),'m--');
% 
% legend([L1,L2,L3,L4,L5],'Perfect MOR L^2 Error','GPE Snapshot bases MOR L^2 Error','Global bases MOR L^2 Error','GPE bases MOR L^2 Error','GPE Snapshot L^2 Error')
% title(sprintf('L^2 error accumulation'))
% hold off


% ----L2 Error accumulation plot Sum of all cases------------------------------------------------------------
% SSE_dxdn_U_OrigSnapS=sum(sum((Y_U_OrigSnapS_Rec-Y_Rec(:,:,Test_StartIndex:Test_StartIndex+Num_Test-1)).^2,3)); %Square sum error; integral on dx & dn(cases)
% SSE_dxdn_U_GlobalSnapS=sum(sum((Y_U_GlobalSnapS_Rec-Y_Rec(:,:,Test_StartIndex:Test_StartIndex+Num_Test-1)).^2,3)); %Square sum error; integral on dx & dn(cases)
% SSE_dxdn_U_GPESnapS=sum(sum((Y_U_GPESnapS_Rec-Y_Rec(:,:,Test_StartIndex:Test_StartIndex+Num_Test-1)).^2,3)); %Square sum error; integral on dx & dn(cases)
% SSE_dxdn_U_GPE=sum(sum((Y_U_GPE_Rec-Y_Rec(:,:,Test_StartIndex:Test_StartIndex+Num_Test-1)).^2,3)); %Square sum error; integral on dx & dn(cases)
% 
% i = 1:Paras.t_n+1;
% x=(i-1)*(Paras.t_end/Paras.t_n);
% figure
% 
% L1=plot(x,cumsum((SSE_dxdn_U_OrigSnapS)),'b--');
% hold on
% L2=plot(x,cumsum((SSE_dxdn_U_GPESnapS)),'r-.');
% L3=plot(x,cumsum((SSE_dxdn_U_GlobalSnapS)),'g--');
% L4=plot(x,cumsum((SSE_dxdn_U_GPE)),'y--');
% % L5=plot(x,cumsum(sum(SE_Y_GPESnapS(:,1:Paras.t_n+1))),'m--');
% 
% legend([L1,L2,L3,L4],'ROM ferfect base','ROM GPE Snapshot base','ROM GPE global base','ROM GPE base')
% % legend([L1,L2,L3,L4,L5],'Perfect MOR L^2 Error','GPE Snapshot bases MOR L^2 Error','Global bases MOR L^2 Error','GPE bases MOR L^2 Error','GPE Snapshot L^2 Error')
% title(sprintf('L^2 error accumulation. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,n_Ubases,dim_new))
% hold off
% set(gca,'yscale','log');
   
% % ----Animation plot------------------------------------------------------
% x = 0:1/Paras.n:1;  % coordinate sequence
% Y_Maxi=max([Y_Full(:);Y_U_OrigSnapS(:);Y_U_GPESnapS(:)]);
% Y_Mini=min([Y_Full(:);Y_U_OrigSnapS(:);Y_U_GPESnapS(:)]);
% for i = 1:Paras.t_n+1
%     figure(5)   
% %     title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))
% %     subplot(2,1,1)
% 
%     L1=plot(x,Y_Full(:,i),'k-');
%     hold on
%     L2=plot(x,Y_U_OrigSnapS(:,i),'b--');
%     L3=plot(x,Y_U_GPESnapS(:,i),'r-.');
%     L4=plot(x,Y_U_GlobalSnapS(:,i),'g--');
%     L5=plot(x,Y_U_GPE(:,i),'y--');
%     
%     if mod(i,Paras.t_n/Paras.n_snap)==1;
% %         index=fix(i/(Paras.t_n/Paras.n_snap));
%         L6=plot(x,Y_GPESnapS(:,fix(i/(Paras.t_n/Paras.n_snap))+1),'m--');
%     end
%         
% %     axis([0,1,Y_Mini,Y_Maxi]);
% %     legend([L1,L2,L3,L4,L5,L6],'Full FEM Solution','Perfect MOR FEM Solution','GPE Snapshot bases MOR FEM Solution','Global bases MOR FEM Solution','GPE bases MOR FEM Solution','GPE prediction')
% %     title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n)))) 
% %     hold off
%     
%     
% % %     figure(5)
% %     subplot(2,1,2)
% %     L1=plot(x(1:i),sum(Y_Full(:,1:i)-Y_Full(:,1:i)),'k-');
% %     hold on
% %     L2=plot(x(1:i),sum(Y_U_OrigSnapS(:,1:i)-Y_Full(:,1:i)),'b--');
% %     L3=plot(x(1:i),sum(Y_U_GPESnapS(:,1:i)-Y_Full(:,1:i)),'r-.');
% %     L4=plot(x(1:i),sum(Y_U_GlobalSnapS(:,1:i)-Y_Full(:,1:i)),'g--');
% %     L5=plot(x(1:i),sum(Y_U_GPE(:,1:i)-Y_Full(:,1:i)),'y--');
%     
%     axis([0,1,Y_Mini,Y_Maxi]);
%     legend([L1,L2,L3,L4,L5],'Full FEM Solution','Perfect MOR FEM Solution','GPE Snapshot bases MOR FEM Solution','Global bases MOR FEM Solution','GPE bases MOR FEM Solution')
%     title(sprintf('L^2 Error Acumulation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n)))) 
%     hold off
% 
%     F(i) = getframe;    
% end    
    





