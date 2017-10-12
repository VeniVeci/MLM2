%% MC_HeatDC_Read
% Read save Monte Carlo process result
%
% Modifications:
% WeiX, 11-4-2016, Create.


clear

% Parameter
%Filter parameter
alpha=2;

%Sensor 
iSensor=5;



load('HeatDC_MC_data13.mat')
[num_sample,~,~]=size(Y_Rec);

% Sum_Y_Rec=sum(Y_Rec,2);
Sum_Y_Rec=trapz(Y_Rec,2)*Paras.dt;

Sum_Y_Rec=Sum_Y_Rec(:,iSensor);


Sum_Y_Rec=reshape(Sum_Y_Rec,num_sample,[]);

% Filter out failure point
SS_Y_Rec=sum(Sum_Y_Rec,2);
med=median(SS_Y_Rec);

index=find(SS_Y_Rec<alpha*med & SS_Y_Rec>=0);

Sensor_HDM=Sum_Y_Rec(index,:);

Time_HDM=sum(Time_Rec);

load('HeatDC_MC_data132.mat')
% ilter out failure point 
% Sum_Y_UGlobal_Rec=sum(Y_UGlobal_Rec,2);
Sum_Y_UGlobal_Rec=trapz(Y_UGlobal_Rec,2)*Paras.dt;

Sum_Y_UGlobal_Rec=reshape(Sum_Y_UGlobal_Rec,num_sample,[]);

%Filter out failure point
SS_Y_UGlobal_Rec=sum(Sum_Y_UGlobal_Rec,2);
med=median(SS_Y_UGlobal_Rec);
index=find(SS_Y_UGlobal_Rec<alpha*med & SS_Y_Rec>=0);

Sensor_UGlobal=Sum_Y_UGlobal_Rec(index,:);


% ilter out failure point
% Sum_Y_GPESnapS_Rec=sum(Y_GPESnapS_Rec,2);
Sum_Y_GPESnapS_Rec=trapz(Y_GPESnapS_Rec,2)*Paras.dt;

Sum_Y_GPESnapS_Rec=reshape(Sum_Y_GPESnapS_Rec,num_sample,[]);

%Filter out failure point
SS_Y_GPESnapS_Rec=sum(Sum_Y_GPESnapS_Rec,2);
med=median(SS_Y_GPESnapS_Rec);
index=find(SS_Y_GPESnapS_Rec<alpha*med & SS_Y_GPESnapS_Rec>=0);

Sensor_UGPEss=Sum_Y_GPESnapS_Rec(index,:);


% Xrange=[1 8]
%% MEAN ,standard deviation and skewness
mean_Sensor_HDM=mean(Sensor_HDM);
mean_Sensor_UGlobal=mean(Sensor_UGlobal);
mean_Sensor_UGPEss=mean(Sensor_UGPEss);

std_Sensor_HDM=std(Sensor_HDM);
std_Sensor_UGlobal=std(Sensor_UGlobal);
std_Sensor_UGPEss=std(Sensor_UGPEss);

skewness_Sensor_HDM=skewness(Sensor_HDM);
skewness_Sensor_UGlobal=skewness(Sensor_UGlobal);
skewness_Sensor_UGPEss=skewness(Sensor_UGPEss);

Time_UGlobal=sum(Time_UGlobal_Rec);
Time_UGPEss=sum(Time_GPESnapS_Rec);

i=5;
report(1,:)=[mean_Sensor_HDM(i),std_Sensor_HDM(i),skewness_Sensor_HDM(i),Time_HDM/3600];
report(2,:)=[mean_Sensor_UGPEss(i),std_Sensor_UGPEss(i),skewness_Sensor_UGPEss(i),Time_UGPEss/3600];
report(3,:)=[mean_Sensor_UGlobal(i),std_Sensor_UGlobal(i),skewness_Sensor_UGlobal(i),Time_UGlobal/3600];

close all
% i=5;
% figure
% histogram(Sensor_HDM(:,i),70)
% title('HDM')
% % xlim(Xrange)
% 
% figure
% histogram(Sensor_UGlobal(:,i),50)
% title('ROM U by global')
% % xlim(Xrange)
% 
% figure
% histogram(Sensor_UGPEss(:,i),50)
% title(sprintf('ROM U by GPRSS. %s, Num_{Trian}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DimNew}=%i', Para_Dr.DrMethod,Para_Dr.Num_Train,Para_Dr.Num_Snapshot,Para_Dr.n_Ubases,Para_Dr.dim_new))
% % xlim(Xrange)




%% Generate EPS figure
i=5;
a=0.002;
b=0.02
num_bar=60;

data1=Sensor_HDM(find(Sensor_HDM(:,i)>a & Sensor_HDM(:,i)<b));
data2=Sensor_UGlobal(find(Sensor_UGlobal(:,i)>a & Sensor_UGlobal(:,i)<b));
data3=Sensor_UGPEss(find(Sensor_UGPEss(:,i)>a & Sensor_UGPEss(:,i)<b));

figure
histogram(data1,num_bar)
title('HDM')

figure
histogram(data2,num_bar)
title('ROM U by global')

figure
histogram(data3,num_bar)
title(sprintf('ROM U by GPRSS'));

% title(sprintf('ROM U by GPRSS. %s, Num_{Trian}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DimNew}=%i', Para_Dr.DrMethod,Para_Dr.Num_Train,Para_Dr.Num_Snapshot,Para_Dr.n_Ubases,Para_Dr.dim_new))


%% Plot Formating
ax = gca; 
ax.FontSize=30;
ax.Position=[0.1,0.13,0.8,0.8];

xlabel({'Accumulated concentration'},'FontSize',30,'FontWeight','bold');
ylabel({'Frequency'},'FontSize',30,'FontWeight','bold');

set(ax, 'xlim', [0.0008 0.0012])
set(ax, 'ylim', [0 500])

ax.XTick = [0.0008:0.0002:0.0012];
ax.YTick = [0:250:500];

grid on
%Adjust windows size 
h = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 700 700];   %Set position and size 
% set(h,'EdgeColor','none')

figure

%% Generate table 

input.data =report;
% Set column labels (use empty string for no label):
input.tableColLabels = {'Mean','Standard deviation','Skewness','Duration'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'FOM','kGPES ROM','Global basis ROM'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.6f',3,'%.2f',1}; % three digits precision for first three columns, one digit for the last
% input.dataFormat = {'%.3f',3}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 1;

% Uses booktabs basic formating rules ('1' = using booktabs, '0' = not
% using booktabs. Note that in order to compile the generated latex output, you have to use the booktabs package
% within your latex document (include it using \usepackage{booktabs}). Also, using the booktabs option cancels the
% borders options.
input.booktabs = 1;


% LaTex table caption:
input.tableCaption = 'MyTableCaption';

% LaTex table label:
input.tableLabel = 'MyTableLabel';

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% call latexTable:
latex = latexTable(input);





