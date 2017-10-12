function [X,U,V,P,options] = Exp2UVPv(Name_field)
% Exp dataset to U V P dataset (vector form)
% Transfer 'experiment' data to UVP 3 filds in long vector matrix
%  U(k,:)=(RecU(:,:,k)(:))  (in a long vector)

load(Name_field);
[nU,dU]=size(RecU(:,:,1));
[nV,dV]=size(RecV(:,:,1));
[nP,dP]=size(RecP(:,:,1));

% Initilization
% clearvars U V P
U=zeros(point,(nU+2)*(dU+2));
V=zeros(point,(nU+2)*(dU+2));
P=zeros(point,(nP)*(dP));


for i = 1:point
    Re=X(i,1);
    u_lid=X(i,2);
%      DisplayUV(RecU(:,:,i),RecV(:,:,i),RecP(:,:,i),Re,tf,dt,lx,ly,nx,ny,u_lid )

    TmepU=zeros(nU+2,dU+2);
    TmepU(:,end)=ones(nU+2,1)*u_lid;
    TmepU(2:end-1,2:end-1)=RecU(:,:,i);
    U(i,:)=(TmepU(:))';
    
    TmepV=zeros(nV+2,dV+2);
%     TmepV(:,end)=ones(nV+2,1);
    TmepV(2:end-1,2:end-1)=RecV(:,:,i);
    V(i,:)=(TmepV(:))';

    TmepP=RecP(:,:,i);
    P(i,:)=(TmepP(:))';
     
end

options.nU=nU+2;
options.dU=dU+2;
options.nV=nV+2;
options.dV=dV+2;
options.nP=nP;
options.dP=dP;

options.dt=dt;
options.lx=lx;
options.ly=ly;
options.nx=nx;
options.ny=ny;

options.tf=tf;
options.dt=dt;
end

