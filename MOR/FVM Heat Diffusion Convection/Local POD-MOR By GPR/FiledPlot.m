function FiledPlot (Paras,T1,T2,T3)


     PosX=Paras.Lx/Paras.Num_node_x/2:Paras.Lx/Paras.Num_node_x:Paras.Lx-Paras.Lx/Paras.Num_node_x/2;
     PosY=Paras.Ly/Paras.Num_node_y/2:Paras.Ly/Paras.Num_node_y:Paras.Ly-Paras.Ly/Paras.Num_node_y/2;
    %  PosX=PosX';
    %  PosY=PosY';
    [CooX,CooY] = meshgrid(PosX,PosY);


     [~,frame]=size(T1);
     
     
if nargin == 2,     
    figure()
    for i=1:frame
        t=(i-1)*Paras.dt;
        T1_i=reshape(T1(:,i),Paras.Num_node_y,Paras.Num_node_x);

        subplot(2,1,1),contour(CooX,CooY,T1_i),
        title(sprintf('Temperature @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        subplot(2,1,2),pcolor(CooX,CooY,T1_i),shading interp,
        title(sprintf('Temperature @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        getframe;
    end
    
elseif nargin == 3
    figure()
    for i=1:frame
        t=(i-1)*Paras.dt;
        T1_i=reshape(T1(:,i),Paras.Num_node_y,Paras.Num_node_x);
        T2_i=reshape(T2(:,i),Paras.Num_node_y,Paras.Num_node_x);
        
        subplot(2,2,1),contour(CooX,CooY,T1_i),
        title(sprintf('T1 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        subplot(2,2,2),pcolor(CooX,CooY,T1_i),shading interp,
        title(sprintf('T1 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        
        subplot(2,2,3),contour(CooX,CooY,T2_i),
        title(sprintf('T2 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        subplot(2,2,4),pcolor(CooX,CooY,T2_i),shading interp,
        title(sprintf('T2 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        getframe;
    end
    
 elseif nargin == 4
    figure()
    for i=1:frame
        t=(i-1)*Paras.dt;
        T1_i=reshape(T1(:,i),Paras.Num_node_y,Paras.Num_node_x);
        T2_i=reshape(T2(:,i),Paras.Num_node_y,Paras.Num_node_x);
        T3_i=reshape(T3(:,i),Paras.Num_node_y,Paras.Num_node_x);
        
        subplot(2,3,1),contour(CooX,CooY,T1_i),
        title(sprintf('T1 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        subplot(2,3,4),pcolor(CooX,CooY,T1_i),shading interp,
        title(sprintf('T1 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        
        subplot(2,3,2),contour(CooX,CooY,T2_i),
        title(sprintf('T2 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        subplot(2,3,5),pcolor(CooX,CooY,T2_i),shading interp,
        title(sprintf('T2 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        
        subplot(2,3,3),contour(CooX,CooY,T3_i),
        title(sprintf('T3 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        subplot(2,3,6),pcolor(CooX,CooY,T3_i),shading interp,
        title(sprintf('T3 @ t=%0.3f',t)),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
        getframe;
    end   
    
    
    
else
    err('not enough input')
    
end
    
    
    
    
    
    
    
    
    
    
end
