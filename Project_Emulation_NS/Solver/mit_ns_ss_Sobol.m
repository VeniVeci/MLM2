
clear 

point=500;

Re_downlimit=500;
Re_uplimit=2500;

Ulid_downlimit=0.01;
Ulid_uplimit=10;

X=i4_sobol_generate ( 2, point, 0 );
X=X';
X(:,1)=X(:,1)*(Re_uplimit-Re_downlimit)+ones(point,1)*Re_downlimit;
X(:,2)=X(:,2)*(Ulid_uplimit-Ulid_downlimit)+ones(point,1)*Ulid_downlimit;



% Re = 10e2;     % Reynolds number
dt = 1e-3;     % time step
tf = 200e-0;    % final time
lx = 1;       % width of box
ly = 1;       % height of box
nx = 100;      % number of x-gridpoints
ny = 100;      % number of y-gridpoints
nsteps = 1000;  % number of steps with graphic output

% u_lid=1;        % speed of lid

SSlimit=0.0005;

% [U,V,P]=  mit_ns_ss(Re,dt,tf,lx,ly,nx,ny,nsteps,u_lid,SSlimit);
% 
% DisplayUV( U,V,P,Re,tf,dt,lx,ly,nx,ny,u_lid )


for i=1:point
    Re=X(i,1);
    u_lid=X(i,2);
    
    
    [U,V,P]=  mit_ns_ss(Re,dt,tf,lx,ly,nx,ny,nsteps,u_lid,SSlimit);
    
    RecU(:,:,i)=U;
    RecV(:,:,i)=V;
    RecP(:,:,i)=P;
    
end

% save('exp2');
clock










