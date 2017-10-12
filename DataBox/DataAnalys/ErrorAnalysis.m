function [l2ErrorStatistic] = ErrorAnalysis(y,yTurth)
%% Compare y and yTurth to indicate error
% y is (dimension*point number) matrix. yTurth must be same size.


%
ADJACENT_WEIGHT=1.5;

%% Analysis for representative case 
sqrError=(y-yTurth).^2;
l2Error=sqrt(sum(sqrError,1));

[l2ErrorStatistic] = boxSort(l2Error);









%% BACK UP CODE
% errorStatistic.mean     = mean(l2Error);
% errorStatistic.min      = min(l2Error);
% errorStatistic.max      = max(l2Error);
% errorStatistic.quantile = quantile(l2Error,[.25 .5 .75]);
% 
% bodyLength=quantile(l2Error,[.75])-quantile(l2Error,[.25]);
% 
% errorStatistic.upAdjacent = quantile(l2Error,[.75])+ ADJACENT_WEIGHT*bodyLength;
% errorStatistic.lowAdjacent= quantile(l2Error,[.25])- ADJACENT_WEIGHT*bodyLength;



end

