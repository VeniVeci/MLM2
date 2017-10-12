function [Y,model] = DataPP(X)
%Data Pre-process
%This function make multivariate data rescaled from -1 to 1 for all
%properties. Resacled data could be recovered.
%
% Modifications:
% WeiX, Jan-2nd-2016, first edition 

%%
[num,dim]=size(X);

model.name='linear rescale';
model.mean=mean(X,1);
model.maxi=max(X);
model.mini=min(X);
model.range=model.maxi-model.mini;
% model.data=X;

% key=median(X);

Y=(X-repmat(model.mean,[num,1]))./repmat(model.range,[num,1]);


end

