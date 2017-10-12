%% mit_ns_ss_v2_Demo


clear

Re = 10e2;     % Reynolds number
dt = 1e-3;     % time step
tf = 80e-0;    % final time
lx = 1;       % width of box
ly = 1;       % height of box
nx = 20;      % number of x-gridpoints
ny = 20;      % number of y-gridpoints
nsteps = 1000;  % number of steps with graphic output

SSlimit=0.0005;
SSlimit=2e-4;

%Constanct B.C.
BC.uN=5; BC.vN=0;
BC.uS=1; BC.vS=0;
BC.uW=0.5; BC.vW=3;
BC.uE=0; BC.vE=1;


% Parametric B.C.
%     x = linspace(0,lx,nx+1); hx = lx/nx;    %temp
%     y = linspace(0,ly,ny+1); hy = ly/ny;
% 
%     % y=a*sin(4*pi.*x).*exp(-b*x);
%     Na=2;    Nb=3;    Nc=1;
%     Sa=2;    Sb=1;    Sc=1;
%     Wa=2;    Wb=1;    Wc=1;
%     Ea=2;    Eb=1;    Ec=1;
% 
% 
%     BC.uN=Na*sin(Nc*pi.*x).*exp(-Nb*x);     BC.vN=0;
%     BC.uS=Sa*sin(Sc*pi.*x).*exp(-Sb*x);     BC.vS=0;
%     BC.uW=0;                                BC.vW=Wa*sin(Wc*pi.*y).*exp(-Wb*y);
%     BC.uE=0;                                BC.vE=Ea*sin(Ec*pi.*y).*exp(-Eb*y);


[U,V,P]=  mit_ns_ss_v2(Re,dt,tf,lx,ly,nx,ny,nsteps,SSlimit,BC);

% figure 
% DisplayUV_v2( U,V,P,Re,tf,dt,lx,ly,nx,ny,BC )

