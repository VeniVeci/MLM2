%% ProjetMOR_Burger1D_FEM_Part2_TestDrGPEMORx3_PODn_DEIM_MCV2_Read
%  
%  Read data saved by ProjetMOR_Burger1D_FEM_Part2_TestDrGPEMORx3_PODn_DEIM_MCV2
%  version 2: using different dataset formate. 
% 
% Modifications:
% 10-Sep-2016, WeiX, first edition 

clear




%% ----------------Load dataset-------------------------------------------
% load('ExpDataV2_27.mat') 
load('Bur_MOR_kGPE_DEIM_Train180Test300SS200DEIMSS200U5to60UDEIM30v73.mat') 
% load('Bur_MOR_kGPE_DEIM_Train180Test300SS200DEIMSS200U5to60UDEIM30.mat') 

h = waitbar(0,'process data');
for i =1:Num_Test
        
        for j=5:5:n_Ubases
            
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


figure
boxplot(RE_U_GPESnapS(:,5:5:n_Ubases))
xlabel('POD basis dimension');
ylabel('Relative error');
title('(c)');

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=0;
dB=-4;

ax.YLim=[10^dB,10^uB];
ytick=logspace(dB,uB,uB-dB+1);

ax.YTick=ytick(1:2:end);
% ax.YTickLabel=ytick(1:2:end);

% ax.XTick=5:10:n_Ubases-5;
ax.XTick=1:2:13;
ax.XTickLabel=5:10:n_Ubases-5;


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








