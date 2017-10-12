function [ yBoxStatistic ] = boxSort(y)
%% boxplot sort 

%% setting 
%constant
ADJACENT_WEIGHT=1.5;
TOP_INDEX_NUMBER=3;

%% Analysis for representative case 
[nRow,mColumn] = size(y);
yBoxStatistic.mean     = mean(y);
yBoxStatistic.min      = min(y);
yBoxStatistic.max      = max(y);
yBoxStatistic.quantile = quantile(y,[.25; .5; .75]);

bodyLength=quantile(y,[.75])-quantile(y,[.25]);

yBoxStatistic.upAdjacent = quantile(y,[.75])+ ADJACENT_WEIGHT*bodyLength;
yBoxStatistic.lowAdjacent= quantile(y,[.25])- ADJACENT_WEIGHT*bodyLength;


yBoxStatistic.topNearMeanIndex=sortByCompareValue(y,yBoxStatistic.mean)
yBoxStatistic.topNearMeanIndex=yBoxStatistic.topNearMeanIndex(1:TOP_INDEX_NUMBER,:);
% yBoxStatistic.topNearMeanValue=y(yBoxStatistic.topNearMeanIndex);

for i=1:mColumn
    yBoxStatistic.topNearMeanValue(:,i)=y(yBoxStatistic.topNearMeanIndex(:,i),i);
end

yBoxStatistic.topNearUpQuantileIndex=sortByCompareValue(y,quantile(y,[.75]));
yBoxStatistic.topNearUpQuantileIndex=yBoxStatistic.topNearUpQuantileIndex(1:TOP_INDEX_NUMBER,:);
% yBoxStatistic.topNearUpQuantileValue=y(yBoxStatistic.topNearUpQuantileIndex);

for i=1:mColumn
    yBoxStatistic.topNearUpQuantileValue(:,i)=y(yBoxStatistic.topNearUpQuantileIndex(:,i),i);
end

yBoxStatistic.topNearUpAdjacentIndex=sortByCompareValue(y,yBoxStatistic.upAdjacent);
yBoxStatistic.topNearUpAdjacentIndex=yBoxStatistic.topNearUpAdjacentIndex(1:TOP_INDEX_NUMBER,:);
% yBoxStatistic.topNearUpAdjacentValue=y(yBoxStatistic.topNearUpAdjacentIndex);

for i=1:mColumn
    yBoxStatistic.topNearUpAdjacentValue(:,i)=y(yBoxStatistic.topNearUpAdjacentIndex(:,i),i);
end


%% BACK UP CODE
% CompareValue=yBoxStatistic.mean;
% [~,index]= sort(abs(y-CompareValue));
% errorStatistic.topNearMeanIndex = index(1:TOP_INDEX_NUMBER);
% errorStatistic.topNearMeanValue = y(index(1:TOP_INDEX_NUMBER));
% 
% CompareValue=yBoxStatistic.mean;
% [~,index]= sort(abs(y-CompareValue));
% errorStatistic.topNearMeanIndex = index(1:TOP_INDEX_NUMBER);
% errorStatistic.topNearMeanValue = y(index(1:TOP_INDEX_NUMBER));


end

function index=sortByCompareValue(y,CompareValue)
    [nRow,mColumn]=size(y);
    [~,index]= sort(abs(y-ones(nRow,1)*CompareValue));

end