function y=Sources(x,a,b)
% sources term function (time independent)
%
% Modifications:
% 30-Sep-2016, WeiX, first edition 
%%
%     y=x;
%     y=sin(pi*x);
%    y=sin(2*pi*x).*exp(x);
%      y=a*sin(b*pi*x);
%      y=1;
%          y=sin(a*pi*x+b);
%     y=a*exp(b*x).*sin(x);
%     y=sin(a*pi*x+b*pi);
%  y=(1-x).*cos(3*pi*a*(x+1)).*exp(-(1+x)*a)+0.5; 
 
    y=a.*exp(b*x);
    
%     y=a*sin(pi.*x).*exp(-b*x)+a*sin(pi.*(x+0.5)).*exp(-b*(x+0.5));
%     y=a*sin(pi.*x).*exp(-b*x)+sin(4*pi.*x);
 
end

    