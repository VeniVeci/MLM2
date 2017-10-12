function [ ux,uy ] = uxy(x,y,Paras)
% Function that define velocitu in a given coordinate
%
% A prior velocity filed is governed by a given equation
%
% Modifications:
% 1-April-2016, WeiX, first edition 

% close all

%Scheme1



%Scheme2

ux=Paras.ampx*(x-Paras.x0)./((x-Paras.x0).^2+(y-Paras.y0).^2);
uy=Paras.ampy*(y-Paras.y0)./((x-Paras.x0).^2+(y-Paras.y0).^2);

end

