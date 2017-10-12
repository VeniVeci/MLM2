%% mit_ns_ss_v2_Sobol


clear 

%% ----------Parameters--------------------

n_run=1000;
dim_run=5;


Re_downlimit=500;
Re_uplimit=1000;

X=i4_sobol_generate ( dim_run, n_run, 0 );
X=X';
X(:,1)=X(:,1)*(Re_uplimit-Re_downlimit)+ones(n_run,1)*Re_downlimit;
% X(:,2)=X(:,2)*(Ulid_uplimit-Ulid_downlimit)+ones(n_run,1)*Ulid_downlimit;

Lid_uplimit=5;
Lid_downlimit=1;
 X(:,2:5)=X(:,2:5)*(Lid_uplimit-Lid_downlimit)+ones(n_run,dim_run-1)*Lid_downlimit;

%% ---------System Parameters---------------
% Re = 10e2;     % Reynolds number
dt = 1e-3;     % time step
tf = 200e-0;    % final time
lx = 1;       % width of box
ly = 1;       % height of box
nx = 50;      % number of x-gridpoints
ny = 50;      % number of y-gridpoints
nsteps = 1000;  % number of steps with graphic output

SSlimit=0.0005;
SSlimit=2e-4;



% [U,V,P]=  mit_ns_ss(Re,dt,tf,lx,ly,nx,ny,nsteps,u_lid,SSlimit); 
% DisplayUV( U,V,P,Re,tf,dt,lx,ly,nx,ny,u_lid )

%% ---Temp-------
x = linspace(0,lx,nx+1); hx = lx/nx;    %temp
y = linspace(0,ly,ny+1); hy = ly/ny;

%% ---------------Main---------
h = waitbar(0,'HDM on the run');
for i=1:n_run
    
    Re=X(i,1);
    
    % 5 parameters
    BC.uN=X(i,2);     BC.vN=0;
    BC.uS=X(i,3);     BC.vS=0;
    BC.uW=0;                                BC.vW=X(i,4);
    BC.uE=0;                                BC.vE=X(i,5);
    
      % 13 parameters
%     Na=X(i,2);    Nb=X(i,3);    Nc=X(i,4);
%     Sa=X(i,5);    Sb=X(i,6);    Sc=X(i,7);
%     Wa=X(i,8);    Wb=X(i,9);    Wc=X(i,10);
%     Ea=X(i,11);    Eb=X(i,12);    Ec=X(i,13);
    
    
    % 9 parameters
%     Na=X(i,2);    Nb=X(i,3);    Nc=1;
%     Sa=X(i,4);    Sb=X(i,5);    Sc=1;
%     Wa=X(i,6);    Wb=X(i,7);    Wc=1;
%     Ea=X(i,8);    Eb=X(i,9);    Ec=1;
    
    % Implement BC
%     BC.uN=Na*sin(Nc*pi.*x).*exp(-Nb*x);     BC.vN=0;
%     BC.uS=Sa*sin(Sc*pi.*x).*exp(-Sb*x);     BC.vS=0;
%     BC.uW=0;                                BC.vW=Wa*sin(Wc*pi.*y).*exp(-Wb*y);
%     BC.uE=0;                                BC.vE=Ea*sin(Ec*pi.*y).*exp(-Eb*y);


    % Implement BC
%     BC.uN=4*Na*sin(Nc*pi.*x).*exp(-Nb*x)+1;     BC.vN=0;
%     BC.uS=4*Sa*sin(Sc*pi.*x).*exp(-Sb*x)+1;     BC.vS=0;
%     BC.uW=0;                                BC.vW=4*Wa*sin(Wc*pi.*y).*exp(-Wb*y)+1; 
%     BC.uE=0;                                BC.vE=4*Ea*sin(Ec*pi.*y).*exp(-Eb*y)+1; 
    
    
    [U,V,P]=  mit_ns_ss_v2(Re,dt,tf,lx,ly,nx,ny,nsteps,SSlimit,BC);
       
    RecU(:,:,i)=U;
    RecV(:,:,i)=V;
    RecP(:,:,i)=P;
    
    waitbar(i/n_run);
    
end
close(h);

% save('exp2');
% save('Exp_v2_2')










