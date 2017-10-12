function [X] = DataPP_Rocv(Y,model)
%Data Pre-process recover
%Correspond to Data Pre-process to recover data. Same data formate.
%
% Modifications:
% WeiX, Jan-2nd-2016, first edition 

if model.name ~= 'linear rescale'
    error('Error input for DataPP_Rocv');
end
[num,dim]=size(Y);

X=Y.*repmat(model.range,[num,1])+repmat(model.mean,[num,1]);


end
