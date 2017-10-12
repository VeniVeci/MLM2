%% ProjetMOR_Burger1D_FEM_Part2_TestGlobalU_PODn_MCV2
%  version 2: using different dataset formate.
% Modifications:
% 10-Sep-2016, WeiX, first edition 

clear

%%---------------Setting parameter---------------------------------------
Num_Train=180;   
Num_Test=300;
Test_StartIndex=201;

Num_Snapshot=40; %100
n_Ubases=30;     %Number of POD bases

% dim_new=6;      % For LPCA more than 10 require. The new dimension of reduced emulatio

%% ----------------Load dataset-------------------------------------------
% load('ExpDataV2_24.mat') %7 is ok for HH, HH fails in 4.
load('ExpDataV2_9.mat') %7 is ok for HH, HH fails in 4.

Index_snapshot=1:(Paras.t_n/Num_Snapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

Y=Y_Rec(:,Index_snapshot,:);
Y=reshape(Y,[],500); %500 for this dataset
Y=Y';
Time_HDM=Time;

X_star=X(Test_StartIndex:Test_StartIndex+Num_Test-1,:);
Y_starorig=Y(Test_StartIndex:Test_StartIndex+Num_Test-1,:);

X=X(1:Num_Train,:);
Y=Y(1:Num_Train,:);

% Paras.Re=X(Index_test,1);
% Paras.u0a=X(Index_test,2);
% Paras.u0b=X(Index_test,3);

%% ---------------HDM----------------------------------------------------
%-Disable As full data is provided.

% tic;
% h = waitbar(0,'Test HDM');
% for i =1:Num_Test
%     Paras.Re=X_star(i,1);
%     Paras.u0a=X_star(i,2);
%     Paras.u0b=X_star(i,3);
%     
%     [Y_Full,~,~]=Burger1D_FEM_DBC_SolverF(Paras);
%     Y_Full_Rec(:,:,i)=Y_Full;
%     waitbar(i/Num_Test);
% end
% close(h);
% Time_HDM=toc;


%% ---------------MOR Bases by perfect snapshot -----------------------------------
% toc0=toc;
% h = waitbar(0,'Test MOR Bases by perfect snapshot');
% for i =1:Num_Test
%     
%     Paras.Re=X_star(i,1);
%     Paras.u0a=X_star(i,2);
%     Paras.u0b=X_star(i,3);
%     Y_orig_snaps=Y_Rec(:,Index_snapshot,Test_StartIndex+i-1);
%     
%     [U_OrigSnapS,~,~]=svd(Y_orig_snaps(2:end-1,:));  % U*S*V'=Rec_X
%     % [U_orig,S,~]=svd(T_orig);
%     U_OrigSnapS=U_OrigSnapS(:,1:n_Ubases);
% %     eigenvalues_OrigSnapS=diag(S);
%     [Y_U_OrigSnapS,~,~]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_OrigSnapS);
%     Y_U_OrigSnapS_Rec(:,:,i)=Y_U_OrigSnapS; 
%     waitbar(i/Num_Test);
% end
% close(h);
% Time_ROM_U_OrigSnapS=toc-toc0;

%% ---------------MOR Golboal Bases -----------------------------------
DrMethod='Global basis';
Y_Global_snaps=reshape(Y',Paras.n+1,[]);
% n_Ubases=5;
% Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point

% [U_GlobalSnapS,~,~]=svd(Y_Global_snaps(2:end-1,:));  % U*S*V'=Rec_X
[U_GlobalSnapS,~,~]=svd(Y_Global_snaps(2:end-1,:),'econ');  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);

% eigenvalues_GlobalSnapS=diag(S);

h = waitbar(0,'Test MOR Bases by Golboal Bases');
for j=1:1:n_Ubases

    U_GlobalSnapSj=U_GlobalSnapS(:,1:j);
    
    for i =1:Num_Test

        Paras.Re=X_star(i,1);
        Paras.u0a=X_star(i,2);
        Paras.u0b=X_star(i,3);
        Y_orig_snaps=reshape(Y_starorig(i,:)',Paras.n+1,[]);

        [Y_U_GlobalSnapS,~,Time_ROM_U_GlobalSnapS]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_GlobalSnapSj);
        
        Y_U_GlobalSnapS_Rec(:,:,i)=Y_U_GlobalSnapS;
        Y_U_GlobalSnapS_Rec2(:,:,i,j)=Y_U_GlobalSnapS;
%         waitbar(i/Num_Test);

    
        SSE_dx_U_GlobalSnapS(i,:)=sum((Y_U_GlobalSnapS_Rec(:,:,i)-Y_Rec(:,:,Test_StartIndex+i-1)).^2,1);   %Square sum error; integral on dx
        SSE_dxdt_U_GlobalSnapS(i,j)=sum(SSE_dx_U_GlobalSnapS(i,:),2);                              %Square sum error; integral on dx & dt
        RE_U_GPESnapS(i,j)=mean(sqrt(SSE_dx_U_GlobalSnapS(i,:) ./ sum(Y_Rec(:,:,Test_StartIndex+i-1).^2,1)));
    end
    waitbar(j/n_Ubases);
 
end
 close(h);

figure
boxplot(SSE_dxdt_U_GlobalSnapS)
title(sprintf('L^2 Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i', DrMethod,Num_Train,Num_Test,Num_Snapshot))
set(gca,'yscale','log');
xlabel('POD basis dimension');
ylabel('Square error');


figure
boxplot(RE_U_GPESnapS)
title(sprintf('Relative Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i', DrMethod,Num_Train,Num_Test,Num_Snapshot))
set(gca,'yscale','log');
xlabel('POD basis dimension');
ylabel('Relative error');




% 
% figure
% boxplot(RE_U_GPESnapS(:,1:1:n_Ubases))
% xlabel('POD basis dimension');
% ylabel('Relative error');
% title('(b)');
% 
% ax = gca; 
% ax.YScale='log';
% %Range for log scale
% ax = gca; 
% uB=0;
% dB=-5;
% 
% ax.YLim=[10^dB,10^uB];
% ytick=logspace(dB,uB,uB-dB+1);
% 
% % ax.YTick=ytick(1:2:end);
% % ax.YTickLabel=ytick(1:2:end);
% 
% % ax.XTick=5:10:n_Ubases-5;
% ax.XTick=5:5:30;
% ax.XTickLabel=5:5:n_Ubases;











% ylabel({'Relative error'},'FontSize',30,'FontName','Times New Roman');
% xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman');
% ax = gca; 
% ax.FontSize=30;

% toc0=toc;
% h = waitbar(0,'Test MOR Bases by Golboal Bases');
% for i =1:Num_Test
%     
%     Paras.Re=X_star(i,1);
%     Paras.u0a=X_star(i,2);
%     Paras.u0b=X_star(i,3);
%     Y_orig_snaps=reshape(Y_starorig(i,:)',Paras.n+1,[]);
%     
%     [Y_U_GlobalSnapS,~,Time_ROM_U_GlobalSnapS]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_GlobalSnapS);
% 
%     Y_U_GlobalSnapS_Rec(:,:,i)=Y_U_GlobalSnapS;
%     waitbar(i/Num_Test);
% end
% close(h);
% Time_ROM_U_GlobalSnapS=toc-toc0;

%  FiledPlot (Paras,T_orig,T_U_orig)


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
    





