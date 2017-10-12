clear

Re = 12e2;     % Reynolds number
dt = 1e-2;     % time step
tf = 80e-0;    % final time
lx = 1;       % width of box
ly = 1;       % height of box
nx = 100;      % number of x-gridpoints
ny = 100;      % number of y-gridpoints
nsteps = 100;  % number of steps with graphic output

u_lid=1;        % speed of lid

SSlimit=0.0005;


% [U,V,P]=  mit_ns_ss(Re,dt,tf,lx,ly,nx,ny,nsteps,u_lid,SSlimit);
[U,V,P]=  mit_ns_ss(Re,dt,tf,lx,ly,nx,ny,nsteps,u_lid,SSlimit);

DisplayUV( U,V,P,Re,tf,dt,lx,ly,nx,ny,u_lid )
