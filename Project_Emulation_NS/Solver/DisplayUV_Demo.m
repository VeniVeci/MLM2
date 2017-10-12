
load('exp1')




% % Re = 10e2;     % Reynolds number
% dt = 1e-3;     % time step
% tf = 200e-0;    % final time
% lx = 1;       % width of box
% ly = 1;       % height of box
% nx = 100;      % number of x-gridpoints
% ny = 100;      % number of y-gridpoints
% nsteps = 1000;  % number of steps with graphic output
% 
% % u_lid=1;        % speed of lid
% 
% SSlimit=0.0005;


for i = 1:500
    i=4;
     Re=X(i,1);
     u_lid=X(i,2);
     DisplayUV(RecU(:,:,i),RecV(:,:,i),RecP(:,:,i),Re,tf,dt,lx,ly,nx,ny,u_lid )
     
end

