function y=U0(x,a,b)
% Initial condition function (time independent)
%
% Modifications:
% 27-May-2015, WeiX, first edition 
% 09-Jan-2017, WeiX, reuse old function. of BC BC %11-Sep-2016
%%
%     y=x;
%     y=sin(pi*x);
%    y=sin(2*pi*x).*exp(x);
%      y=a*sin(b*pi*x);
%    y=a*exp(b*x);

%      y=sin(a*pi*x+b);       %probably latest used
%     y=a*exp(b*x).*sin(x);
%     y=sin(a*pi*x+b*pi);
%     y=a*exp(b*x).*sin(pi*x);    % continous with Dirichlet BC %10-Sep-2016
    y=1*exp(-(a*x+b)).*sin(pi*x);    % continous with Dirichlet BC %11-Sep-2016
%     y=1*exp(-(a*x+b)).*sin(3*pi*x);    % continous with Dirichlet BC %1-Oct-2016 
end

    