function [Result] = CmpYField(Y_star,Y,options)
% Description:
% Compare Y and Y_origi Field vector
%
% Synopsis:
% [Result] = CmpYField(Y,Y_origi)
% [Result] = CmpYField(Y,Y_origi,options)
%
% Input: 
%       Y_star         % [samples X dimension X fold] matrix, data to be compared
%       Y              % [samples X dimension] matrix, standard data 
%       options        % options of output
%
% Output: 
%       Result         % Result
%
%  Modification
%  WeiX, 6-1-2016,  create
%
%% Initialization and Parameters
% [num,dim]=size(X);
if nargin < 2, options = []; end
if ~isfield(options,'chart'), options.chart = 1; end     

%% Main
[num,dim,fold]=size(Y_star);

for i = 1: fold

    SE=(Y_star-Y).^2;                   %Square error of Y^(i)_j
    
    SS_Y(:,i)=sum(Y.^2,2);              %Square sum of Y^(i)
    MS_Y(:,i)=mean(Y.^2,2);             %Mean of square
     
    SSE(:,i)=sum(SE,2);                 %Square sum error of Y^(i)
    [SSE_sr(:,i),SSE_index(:,i)]=sort( SSE(:,i));
    
    RSSE(:,i)=mean(SE,2)./MS_Y(:,i);    %Relative SSE
    [RSSE_sr,RSSE_index(:,i)]=sort( RSSE(:,i));

%     SRSE(:,i)=sum(SE./(Y.^2),2);      %Sum of relative square error SRSE=SSE
%     [SRSE_sr,SRSE_index(:,i)]=sort( SRSE(:,i));

end

%!! for many time 

Result.SSE=SSE;
Result.SSE_index=SSE_index;

Result.RSSE=RSSE;
Result.RSSE_index=RSSE_index;

% Result.SRSE=SRSE;
% Result.SRSE_index=SRSE_index;

%% Plot SSE relative 
% plot only available when options.fold=1

figure
boxplot(SSE)
title('SSE')

figure
plot(SSE_sr)
title('SSE (Rearragened index)')

figure 
plot(Y_star(SSE_index(1),:),'k-');
hold on 
plot(Y(SSE_index(1),:),'r--');
hold off
legend('Y_star','Y')
title('Best case according to SSE')

figure 
plot(Y_star(SSE_index(end),:),'k-');
hold on 
plot(Y(SSE_index(end),:),'r--');
hold off
legend('Y_star','Y')
title('Worst case according to SSE')

figure 
plot(Y_star(SSE_index(round(num/2)),:),'k-');
hold on 
plot(Y(SSE_index(round(num/2)),:),'r--');
hold off
legend('Y_star','Y')
title('Median case according to SSE')



%% Plot RSSE relative
figure
boxplot(RSSE)
title('RSSE')

figure
plot(RSSE_sr)
title('RSSE (Rearragened index)')

figure 
plot(Y_star(RSSE_index(1),:),'k-');
hold on 
plot(Y(RSSE_index(1),:),'r--');
hold off
legend('Y_*','Y')
title('Best case according to RSSE')

figure 
plot(Y_star(RSSE_index(end),:),'k-');
hold on 
plot(Y(RSSE_index(end),:),'r--');
hold off
legend('Y_*','Y')
title('Worst case according to RSSE')

figure 
plot(Y_star(RSSE_index(round(num/2)),:),'k-');
hold on 
plot(Y(RSSE_index(round(num/2)),:),'r--');
hold off
legend('Y_*','Y')
title('Median case according to RSSE')











