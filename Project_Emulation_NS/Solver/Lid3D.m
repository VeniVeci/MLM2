%LDC3D8HW   Lid-driven cavity 
% Finite element solution of the 3D Navier-Stokes equation 
%   using 8-node, 48-DOF linear rectangular velocity 
%   basis for flow in a rectangular 3D lid-driven cavity.
% The rectangular problem domain is defined between Cartesian 
%   coordinates Xmin & Xmax, Ymin & Ymax and Zmin & Zmax.
% Uses symmetry about z=0 to model one-half of cube.
% The computational grid has NumEx elements in the x-direction, 
%   NumEy elements in the y-direction and NumEz elements in the z-direction.
%
% References:
%    Gartling and Reddy, FEM in Heat Transfer and Fluid Dynamics, p196-200.
%    Ding, Shu, Yeo & Xu,  Comp Meth Appl Mech & Engr 195 (2006) 516-533.
%    Jiang, Lin & Provinelli, Comp Meth Appl Mech & Engr 114 (1994) 213-231. 
%    Holdeman, J.T., A velocity-stream function method for three-dimensional incompressible 
%      fluid flow, Comp Meth Appl Mech & Engr, (conditionally accepted July, 2011). 
%
% Calls: 
%   DMat3D8W            - Element diffusion matrix
%   CMat3D8W           - Element convection matrix 
%   regrade             - for graded mesh
%   ilu_gmres_with_EBC  - Equation solver (GMRES wirh ILU preconditioning)
%
% Include:
%   V8cW       - Velocity element for use in convection matrix
%   V8xyzW     - Velocity derivatives for use in diffusion & convection matricies

% Jonas Holdeman - January 2011,  revised July 2011 

clear all;

disp('3D lid-driven cavity, half of cube.');
disp(' Linear 8-node divergence-free elements.');
disp(' ');

% ---------------------------------------------------------
nd = 6; nd2=nd*nd;   % Number DOFs per node. DO NOT CHANGE!
nv = 8; nv2=nv*nv;   % Number of nodes per element. 
% ---------------------------------------------------------
ETstart=clock;

% Set parameters for GMRES solver 
GMRES.Tolerance=1.e-12;
GMRES.MaxIterates=14; 
GMRES.MaxRestarts=4;

% Set mesh bounds 
Xmin = 0;   Xmax = 1;   DX=Xmax-Xmin;
Ymin = 0;   Ymax = 1;   DY=Ymax-Ymin;
Zmin = 0;   Zmax = .5;  DZ=Zmax-Zmin; % full [-1/2.1/2]
% Calculate hydraulic diameter (for rectangular duct) = 4*area/perimeter
Lc=DY*DZ/(DY+DZ);   % Characteristic length 

% Set mesh grading parameters
xgrd = 1; ygrd=1; zgrd=1; % 
%xgrd = .80; ygrd=.80; zgrd=.80; % graded 

% Set number of elements 
NumEx = 14; 
NumEy = 14; 
NumEz = 7;%   % half of mesh
NumEL = NumEx*NumEy*NumEz;

% Calculate number of nodes 
NumNx=NumEx+1;
NumNy=NumEy+1;
NumNz=NumEz+1;

% Calculate maximum number of nodes 
NumNod=NumNx*NumNy*NumNz;
% Calculate maximum number of degrees-of-freedom 
MaxDof=nd*NumNod;

% Set mean flow velocity (x-direction) 
Ulid=1;
% Set assumed Reynolds number
Re=100;
% Calculate viscosity for specified Reynolds number
nu=Ulid*DY/Re/1.0;  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< correct this!!!!!

% factor for under-relaxation starting at iteration RelxStrt 
% Suggest RelxFac=.9 for Re=100, RelxFac=.5 for Re=400 
RelxFac = .9; 
RelxStrt = 1;

% Set number of nonlinear iterations
MaxNLit=5; %36; %

disp(['Number of elements = ' num2str(NumEL) ', Number of nodes = ' num2str(NumNod) ', Max DOF = ' num2str(MaxDof)]);
disp(['Maximum number of nonlinear iterations = ' num2str(MaxNLit) 'Re = ' num2str(Re)]);
pause(1);

% Mesh generation, corner node coordinates (half mesh)
XN = linspace(Xmin,Xmax,NumNx);
YN = linspace(Ymin,Ymax,NumNy);
ZN = linspace(Zmin,Zmax,NumNz);

% Scale coordinates if graded mesh
if xgrd ~= 1 XN=regrade(XN,xgrd,0); end
if ygrd ~= 1 YN=regrade(YN,ygrd,0); end
if zgrd ~= 1 ZN=regrade(ZN,zgrd,0); end

% Allocate storage for fields (visualization only)
a0=zeros(NumNx,NumNy,NumNz);
b0=zeros(NumNx,NumNy,NumNz);
c0=zeros(NumNx,NumNy,NumNz);
u0= zeros(NumNx,NumNy,NumNz);
v0= zeros(NumNx,NumNy,NumNz);
w0= zeros(NumNx,NumNy,NumNz);

%--------------------Begin grid plot-----------------------
% ********************** FIGURE 1 *************************
% Plot the grid 
figure(1);
clf;
orient portrait;
subplot(2,2,1);
hold on;
  for m=[1:NumNx] plot3([ZN(NumNz);ZN(NumNz)],[XN(m);XN(m)],[YN(1);YN(NumNy)],'k'); end % Face Zmax y
  for m=[1:NumNy] plot3([ZN(NumNz),ZN(NumNz)],[XN(1),XN(NumNx)],[YN(m),YN(m)],'k'); end %      Zmax x
  for m=[1:1] plot3([ZN(1),ZN(1)],[XN(1),XN(NumNx)],[YN(m),YN(m)],':k'); end        %      Zmin x
  for m=[1:NumNx] plot3([ZN(1),ZN(NumNz)],[XN(m),XN(m)],[YN(NumNy),YN(NumNy)],'k'); end % Face Ymax z
  for m=[1:NumNz] plot3([ZN(m),ZN(m)],[XN(1),XN(NumNx)],[YN(NumNy),YN(NumNy)],'k'); end %      Ymax x
  for m=[1:1] plot3([ZN(1),ZN(NumNz)],[XN(1),XN(1)],[YN(m),YN(m)],':k'); end            % Face Xmin z
  for m=[1:1] plot3([ZN(m),ZN(m)],[XN(1),XN(1)],[YN(1),YN(NumNy)],':k'); end            %      Xmin y
  for m=[1:NumNz] plot3([ZN(m),ZN(m)],[XN(NumNx),XN(NumNx)],[YN(1),YN(NumNy)],'k'); end %Face Xmax y
  for m=[1:NumNy] plot3([ZN(1),ZN(NumNz)],[XN(NumNx),XN(NumNx)],[YN(m),YN(m)],'k'); end %Face Xmax z
hold off;
xlabel('z');  ylabel('x');  zlabel('y','Rotation',90);
axis([Zmin,Zmax,Xmin,Xmax,Ymin,Ymax]);
view(140,12);
axis image;
title([num2str(NumEx) 'x' num2str(NumEy) 'x' num2str(NumEz) ...
      ' element mesh for lid-driven cavity']);
  
pause(1);
%-------------- End plotting Figure 1 ----------------------

NodNdx=zeros(3,NumNod); % node number -> (nx,ny,nz)
NodLst=zeros(3,NumNod); % node number -> (x,y,z)
nn2nft=zeros(2,NumNod); % node number -> nf & nt
% nt flags nodal data type (3 velocities or vector potentials): 
%   nt=1 -> (3 vector potential, 3 velocities) 
ni=0; nf=0; dnf=6; nt=1;
for nz=1:NumNz
  for ny=1:NumNy  
    for nx=1:NumNx   
      ni=ni+1;  NodNdx(:,ni)=[nx;ny;nz]; 
      NodLst(:,ni)=[XN(nx);YN(ny);ZN(nz)];
      nn2nft(:,ni)=[nf+1;nt]; nf=nf+dnf;
    end  % loop on nx 
  end   % loop on ny 
end   % loop on nz

% Compute node number (nn) and freedom type (ft) for each freedom number (nf)

nf2nnt=zeros(2,MaxDof);  % node & type associated with dof
ndof=0;
for n=1:NumNod
    ndof=ndof+1;  nf2nnt(:,ndof)=[n;1];  % a 
    ndof=ndof+1;  nf2nnt(:,ndof)=[n;2];  % b 
    ndof=ndof+1;  nf2nnt(:,ndof)=[n;3];  % c 
    ndof=ndof+1;  nf2nnt(:,ndof)=[n;4];  % u 
    ndof=ndof+1;  nf2nnt(:,ndof)=[n;5];  % v 
    ndof=ndof+1;  nf2nnt(:,ndof)=[n;6];  % w 
end

NEx=NumEx; NEy=NumEy; NEz=NumEz;
Elcon = zeros(8,NumEL); 
% GENERATE ELEMENT CONNECTIVITY on block (NEx)x(NEy)x(NEz) 
% Elements are generated increasing along the x-axis, then y-axis, then z-axis. 
% Assumes nodes are generated increasing along the x-axis, then y-axis, then z-axis. 
%                                                ____________ 
%  Element nodal order:                         /|          /|     y 
% [[-1,-1,-1],[1,-1,-1],[-1,1,-1],[1,1,-1],    /___________/ |     |
%  [-1,-1,1],[1,-1,1],[-1,1,1],[1,1,1]]        | |         | |     |____ x
%                                              |           |--->   /
%                                              | |_  _  _  | |    /
%                                              | /         | /   z
%                                              |___________|/
%
ne=0;
Lx =NEx+1; Ly=NEy+1;
Lxy=Lx*Ly;
for nz=1:NEz
  for ny=1:NEy
    for nx=1:NEx
      ne=ne+1;
      Elcon(1,ne)=1+(nx-1)+(ny-1)*Lx+(nz-1)*Lxy;
      Elcon(2,ne)=1+(nx)+(ny-1)*Lx+(nz-1)*Lxy;
      Elcon(3,ne)=1+(nx-1)+(ny)*Lx+(nz-1)*Lxy;
      Elcon(4,ne)=1+(nx)+(ny)*Lx+(nz-1)*Lxy;
      Elcon(5,ne)=1+(nx-1)+(ny-1)*Lx+(nz)*Lxy;
      Elcon(6,ne)=1+(nx)+(ny-1)*Lx+(nz)*Lxy;
      Elcon(7,ne)=1+(nx-1)+(ny)*Lx+(nz)*Lxy;
      Elcon(8,ne)=1+(nx)+(ny)*Lx+(nz)*Lxy;
    end
  end
end
if (NumEL>ne) Elcon=Elcon(:,1:ne); NumEL=ne; end  % trim if necessary

%return;

% Begin ESSENTIAL (Dirichlet) boundary conditions 
%MaxEBC = 2*(NEy+NEz)*(5*NEx+4)+2*5*NEy*NEz;
MaxEBC = 6*2*(NEx*NEy+NEy*NEz+NEz*NEx);
EBC.dof=zeros(MaxEBC,1);  % Degree-of-freedom index  
EBC.val=zeros(MaxEBC,1);  % Dof value 
EBC.num=MaxEBC;
EBC.MxF=MaxDof;

% Simulate full cube, 
nc=0; nU=0; 
for nf=1:MaxDof
  ni=nf2nnt(1,nf);   % which node?
  x=NodLst(1,ni); y=NodLst(2,ni); z=NodLst(3,ni);
  
  if (x==Xmin | x==Xmax | y==Ymin) % Sides & bottom
     switch nf2nnt(2,nf);   % which type?
    case 1   % Ax       
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % a
    case 2   % Ay 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % b
    case 3   % Az 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % c
    case 4   % u 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % u
    case 5   % v 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % v
    case 6   % w 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % w
    end   % switch (type)

  elseif (z==Zmax) % Outside
     switch nf2nnt(2,nf);   % which type?
    case 1   % Ax       
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % a
    case 2   % Ay 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % b
    case 3   % Az 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % c
    case 4   % u 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % u
    case 5   % v 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % v
    case 6   % w 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % w
    end   % switch (type)
    
  elseif (y==Ymax)    % Top (+y is up)
    switch nf2nnt(2,nf);
    case 1   % Ax 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % a
    case 2   % Ay 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % b
    case 3   % Az 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % c
    case 4   % u 
        nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=Ulid; % u*
    case 5   % v 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % v
    case 6   % w 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % w
    end   % switch (type)
    
  elseif (z==Zmin)   % Inside
    switch nf2nnt(2,nf);
    case 1   % Ax 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % a
    case 2   % Ay 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % b
    case 6   % w 
      nc=nc+1;     EBC.dof(nc)=nf; EBC.val(nc)=0;     % w
    end
 
  end   % if 
end   % for nf

disp(['nc= ' num2str(nc) ', MaxEBC= ' num2str(MaxEBC)]);

EBC.num=nc;
if (size(EBC.dof,1)>nc)   % Truncate arrays if necessary 
   EBC.dof=EBC.dof(1:nc);
   EBC.val=EBC.val(1:nc);
end 
% End ESSENTIAL (Dirichlet) boundary conditions 

% partion out essential (Dirichlet) dofs
p_vec = [1:EBC.MxF]';         % List of all dofs
EBC.p_vec_undo = zeros(1,EBC.MxF);
% form a list of non-diri dofs
EBC.ndro = p_vec(~ismember(p_vec, EBC.dof));	% list of non-diri dofs
% calculate p_vec_undo to restore Q to the original dof ordering
EBC.p_vec_undo([EBC.ndro;EBC.dof]) = [1:EBC.MxF]; %p_vec';

Q=zeros(MaxDof,1); % Allocate space for solution (dof) vector

% Initialize dofs to boundary conditions
for k=1:EBC.num
  Q(EBC.dof(k))=EBC.val(k);
end

% Arrays for convergence norm info
MxNL=max(1,MaxNLit);
naa=zeros(1,MxNL);     % Arrays for convergence info
nv0=zeros(1,MxNL);

Dmat = spalloc(MaxDof,MaxDof,80*MaxDof);   % to save the diffusion matrix
Vdof=zeros(6,8);
Xe=zeros(3,8);      % coordinates of element corners 

ItType=0;
NLitr=0; 

% >>>>>>>>>>> BEGIN NONLINEAR ITERATION >>>>>>>>>>>>>>>>>>>>>>  
while (NLitr<MaxNLit), NLitr=NLitr+1;   
      
tclock=clock;   % Start assembly time <<<<<<<<<
% Generate and assemble element matrices
Mat=spalloc(MaxDof,MaxDof,80*MaxDof);
RHS=spalloc(MaxDof,1,MaxDof);
Emat=zeros(48*48,1);         % Values (nv*nv)*(nd*nd) 
DEmat=zeros(48*48,1);    
Rndx=zeros(48*48,1);     
Cndx=zeros(48*48,1);       

% BEGIN GLOBAL MATRIX ASSEMBLY
for ne=1:NumEL   
  for k=1:8
    ki=NodNdx(:,Elcon(k,ne));
    Xe(:,k)=[XN(ki(1));YN(ki(2));ZN(ki(3))]; 
  end  % loop (corner nodes)
  
  if NLitr == 1     
%     Fluid element diffusion matrix, save on first iteration    
    [DEmat,Rndx,Cndx] = DMat3D8W(Xe,Elcon(:,ne),nn2nft);
    Dmat=Dmat+sparse(Rndx,Cndx,DEmat,MaxDof,MaxDof);  % Global diffusion matrix 
  end 

  if (NLitr>1) 
%     Fluid element convection matrix, first iteration uses Stokes equation. 
% Get vector potentials and velocities 
    for k=1:8    % Loop over local nodes of element
      ni=Elcon(k,ne);      % node 
      nf=nn2nft(1,ni);     % dof number/index 
      Vdof(1:6,k)=Q(nf:nf+5);  % u,v,w 
    end    % loop (element nodes)
 
    [Emat,Rndx,Cndx] = CMat3D8W(Xe,Elcon(:,ne),nn2nft,Vdof);    % Convection term for simple iteration
    Mat=Mat+sparse(Rndx,Cndx,-Emat,MaxDof,MaxDof);  % Global convection assembly 
  end  % if (NLitr>1 )

end   % loop ne over elements 
% END GLOBAL MATRIX ASSEMBLY

Mat = Mat - nu*Dmat;    % Add in cached/saved global diffusion matrix 

disp(['(' num2str(NLitr) ')  Matrix assembly complete, elapsed time = '...
    num2str(etime(clock,tclock)) ' sec.  Start solution.']);  % Assembly time <<<<<<<<<<<
pause(1);

Q0 = Q; 

% Solve system
tclock=clock;  % Start solution time  <<<<<<<<<<<<<<

RHSr=RHS(EBC.ndro)-Mat(EBC.ndro,EBC.dof)*EBC.val;
Matr=Mat(EBC.ndro,EBC.ndro);

% condest(Matr)

Qs=Q(EBC.ndro);

Qr=ilu_gmres_with_EBC(Matr,RHSr,[],GMRES,Qs);

Q=[Qr;EBC.val];         % Augment active dofs with esential (Dirichlet) dofs
Q=Q(EBC.p_vec_undo);       % Restore natural order
   
stime=etime(clock,tclock); % Solution time <<<<<<<<<<<<<<

% ****** RELAXATION FACTOR ***************************
if(NLitr>RelxStrt) Q=RelxFac*Q+(1-RelxFac)*Q0; end
% ****************************************************

% Compute change and copy dofs to field arrays
dsqa=0; dsqv=0;
for k=1:MaxDof
  nt=nf2nnt(2,k);
  if     (1<=nt & nt<=3) dsqa=dsqa+(Q(k)-Q0(k))^2;
  elseif (4<=nt & nt<=6) dsqv=dsqv+(Q(k)-Q0(k))^2;
  end  % if  (types)
end  % loop (dofs)

% Copy solution to field arrays
for k=1:MaxDof
  Ni=nf2nnt(1,k);   % which node?
  nx=NodNdx(1,Ni);  ny=NodNdx(2,Ni);  nz=NodNdx(3,Ni);
  switch nf2nnt(2,k)
  case 1
    a0(nx,ny,nz)=Q(k);
  case 2
    b0(nx,ny,nz)=Q(k);
  case 3
    c0(nx,ny,nz)=Q(k);
  case 4
    u0(nx,ny,nz)=Q(k);
  case 5
    v0(nx,ny,nz)=Q(k);
  case 6
    w0(nx,ny,nz)=Q(k);
  end  % switch (type)
end  % loop (dofs)
naa(NLitr)=sqrt(dsqa); 
nv0(NLitr)=sqrt(dsqv); 

disp(['(' num2str(NLitr) ')  Solution time for linear system = ' num2str(etime(clock,tclock)) ...
    ' sec,' ' dV = ' num2str(nv0(NLitr)) ', dA = ' num2str(naa(NLitr))]); % Solution time <<<<<<<
pause(1);

%---------- Begin plot of intermediate results ----------
% ********************** FIGURE 1 *************************
figure(1);

subplot(2,2,3);
nxc=fix((NumNx+1)/2); nzc=1; % nzc=fix((NumNz+1)/2); 
Uc=u0(nxc,:,nzc);
axis([-.3*Ulid,Ulid,Ymin,Ymax]);
plot(Uc,YN,'k+');  % Plot contours (trajectories)
hold on;
plot([-.3,1],[0,1],'w');
hold off;
xlabel('u-velocity');
ylabel('y-coordinate');
title(['Vertical centerline velocity,  Re=' num2str(Re)]);
axis image;
axis square;

subplot(2,2,4);
nyc=fix((NumNy+1)/2); nzc=1; % nzc=fix((NumNz+1)/2); 
Vc=v0(:,nyc,nzc);
plot(XN,Vc,'k+');  % Plot contours (trajectories)
hold on;
plot([0,1],[-.5,.3],'w');
hold off;
axis([Xmin,Xmax,-.5,.3]);
xlabel('x-coordinate');
ylabel('v-velocity');
title(['Horizontal centerline velocity,  Re=' num2str(Re)]);
axis image;
axis square;
 
subplot(2,2,2);
semilogy(1:NLitr,nv0(1:NLitr),'k-+',1:NLitr,naa(1:NLitr),'k-o');  %************
xlabel('Nonlinear iteration number');
ylabel('Nonlinear correction');
title(['Iteration conv.,  Re=' num2str(Re)]);
legend('U','VP');
axis square;

% end figure(1);
%---------- End plot of intermediate results ----------
pause(1);

if (nv0(NLitr)<1e-15) break; end  % Terminate iteration if non-significant 

end;   % <<< END NONLINEAR ITERATION

format short g;
disp('Convergence results by iteration: velocity, vector potential');
disp(['nv0:  ' num2str(nv0(1:NLitr))]); disp(['naa:  ' num2str(naa(1:NLitr))]); 

% ------------------Begin Figure 2 -----------------------

figure(2);
clf;
orient portrait;
orient tall;

subplot(2,2,1);
nc=1;  % fix((NumNz+1)/2);  % plot plane
zp=ZN(nc);
Up=squeeze(u0(:,:,nc));
Vp=squeeze(v0(:,:,nc));
quiver(XN,YN,Up',Vp',2,'k');  % Plot vector field 
hold on;
plot([Xmin,Xmin,Xmax,Xmax,Xmin],[Ymax,Ymin,Ymin,Ymax,Ymax],'k');
hold off;
axis image;
xlabel('x');  ylabel('y'); 
title(['Velocity field in plane z= ' num2str(zp)]);

subplot(2,2,2);
ndy=fix((NumNy+1)/2);
yp=YN(ndy);
Up=squeeze(u0(:,ndy,:));
Wp=squeeze(w0(:,ndy,:));
quiver(ZN,XN,Wp,Up,'k');  % Plot vector field 
hold on;
plot([Zmin,Zmin,Zmax,Zmax,Zmin],[Xmax,Xmin,Xmin,Xmax,Xmax],'k');
hold off;
axis image;
xlabel('z');  ylabel('x'); 
title(['Velocity field in plane y= ' num2str(yp)]);

subplot(2,2,3);
ndx=fix((NumNx+1)/2);
xp=XN(ndx);
Vp=squeeze(v0(ndx,:,:));
Wp=squeeze(w0(ndx,:,:));
quiver(ZN,YN,Wp,Vp,'k');  % Plot vector field 
hold on;
hold off;
axis image;
xlabel('z');  ylabel('y'); 
title(['Velocity field in plane x= ' num2str(xp)]);

% Plot some stream lines 
subplot(2,2,4);
[x0,y0,z0]=meshgrid(XN,YN,ZN);
u1=permute(u0,[2,1,3]);  % Interchange columns for plotting 
v1=permute(v0,[2,1,3]);
w1=permute(w0,[2,1,3]);
sx0=[.7,.25];            % Starting points for stream lines 
sy0=[.84,.5];
sz0=[.02,.38];
 h0=streamline(x0,y0,z0,u1,v1,w1,sx0,sy0,sz0);         
 set(h0, 'Color', 'red');                     
 daspect([1 1 1])                            
 axis tight; box on                          
 axis([Xmin,Xmax,Ymin,Ymax,Zmin,Zmax]);
 xlabel('x');  ylabel('y');  zlabel('z'); 
 camproj perspective;       % perspective or orthographic
 camva('auto')              % Camera viewing angle
 campos([-4 2 4]);          % Camera position 
 camtarget([.5 .5 .25])     % Aiming point 
 camup([0,1,0]);            % Point z-axis up 
 
% ------------------ End Figure 2 ----------------------- 

disp(['Total elapsed time = '...
   num2str(etime(clock,ETstart)/60) ' min']); % Elapsed time from start <<<<<
return;