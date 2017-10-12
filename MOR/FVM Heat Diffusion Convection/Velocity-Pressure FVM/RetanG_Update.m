function [Grid] = RetanG_Update(Grid,x,type)
% retangle grid field update for finite volumn method. 
%    
%
% Modifications:
% 12-May-2014, WeiX, first edition 

switch type
    case 'u'
        u=x;
        [rowu,colu]=size(u);
        [rowGu,colGu]=size(Grid.u);
        
        if rowu==rowGu & colu==colGu
            Grid.u=u;
            %-------interpret glovalu field using linear combination of u--------------
            for x=1:2:Grid.Info.Num_x_full
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalu(y,x)=u(y/2,(x+1)/2);
                end
            end
            for x=1:2:Grid.Info.Num_x_full
                for y = 3:2:Grid.Info.Num_y_full-2
                    Grid.globalu(y,x)=(Grid.globalu(y-1,x)+Grid.globalu(y+1,x))/2;
                end
            end
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalu(y,x)=(Grid.globalu(y,x-1)+Grid.globalu(y,x+1))/2;
                end
            end
        else
            error('size of u doed not match')
        end
    case 'v'
        v=x;
        [rowv,colv]=size(v);
        [rowGv,colGv]=size(Grid.v);   
        
        if rowv==rowGv & colv==colGv
            
            Grid.v=v;
            %-------interpret glovalv field using linear combination of v--------------
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 1:2:Grid.Info.Num_y_full
                    Grid.globalv(y,x)=v((y+1)/2,x/2);
                end
            end
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalv(y,x)=(Grid.globalv(y-1,x)+Grid.globalv(y+1,x))/2;
                end
            end
            for x=3:2:Grid.Info.Num_x_full-2
                for y = 1:2:Grid.Info.Num_y_full
                    Grid.globalv(y,x)=(Grid.globalv(y,x-1)+Grid.globalv(y,x+1))/2;
                end
            end
            
        else
            error('size of v doed not match')
        end

    case 'p'
        p=x;
        [rowp,colp]=size(p);
        [rowGp,colGp]=size(Grid.p);   
        
        if rowp==rowGp & colp==colGp
            
            Grid.p=p;       
            %-------update pressure field (p) in a global grid -----------------------
            for x=2:2:Grid.Info.Num_x_full-1
                for y = 2:2:Grid.Info.Num_y_full-1
                    Grid.globalp(y,x)=p(y/2,x/2);
                end
            end
            
        else
            error('size of p doed not match')
        end
        
    otherwise
        error('No such type in UniGrid');
end

end