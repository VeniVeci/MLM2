% Experiment data ReadMe file
NS_Lid_Experiment1 run on script 'mit_ns_ss_Sobol.m' with condition as follow:
CONDITIONS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% exp1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
point=500;

Re_downlimit=700;
Re_uplimit=1200;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% exp2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%