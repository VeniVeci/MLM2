function [Grid] = RetanG_correction(Grid,p_corrector,alpha_p)
%% retangle grid correction for field using P_corrector
%
% Modifications:
% 12-May-2014, WeiX, first edition 

%%
dx=Grid.dx;
dy=Grid.dy;
Num_x_u=Grid.Num_x_u;
Num_y_u=Grid.Num_y_u;
Num_x_v=Grid.Num_x_v;
Num_y_v=Grid.Num_y_v;

Grid.p=Grid.p+p_corrector*alpha_p;

%check if the p filed is correct

%---method 1------------
% minip=min(Grid.p(:));
% if minip<0
% Grid.p=Grid.p-minip*ones(Grid.Num_x_p,Grid.Num_y_p);
% end

%---method 2------------
% for x= 1: Num_x_p
%     for y= 1: Num_y_p    
%      if p_star(x,y) < 0 p_star(x,y)=0;end
%     end
% end

% ----------u field-------------------------------
for x= 2: Num_x_u-1
    for y= 1: Num_y_u
        p_star_w=p_corrector(x-1,y);
        p_star_e=p_corrector(x,y);
        Grid.u(x,y)=Grid.u(x,y)+(p_star_w-p_star_e)*dy/Grid.a_P_u(x,y);     
    end
end
% ----------v field-------------------------------
for x= 1: Num_x_v
    for y= 2: Num_y_v-1
        
        p_star_s=p_corrector(x,y-1);
        p_star_n=p_corrector(x,y);
        Grid.v(x,y)=Grid.v(x,y)+(p_star_s-p_star_n)*dx/Grid.a_P_v(x,y);     

    end
end



end