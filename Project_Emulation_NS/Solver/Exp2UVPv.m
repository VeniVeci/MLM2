function [X,U,V,P,options] = Exp2UVPv(Name_field)
% Exp dataset to U V P dataset (vector form)
% Transfer 'experiment' data to UVP 3 filds in long vector matrix
%  U(k,:)=(RecU(:,:,k)(:))  (in a long vector)
% RecU is the saved U field.


% Default input parameter
% if nargin <=1
%     number_train=round(point/5);
%     num_test=point-number_train;
% end
    
load(Name_field);
[nxU,nyU]=size(RecU(:,:,1));
[nxV,nyV]=size(RecV(:,:,1));
[nxP,nyP]=size(RecP(:,:,1));

% Initilization
% clearvars U V P
U=zeros(point,(nxU+2)*(nyU+2));
V=zeros(point,(nxU+2)*(nyU+2));
P=zeros(point,(nxP)*(nyP));


for i = 1:point
    Re=X(i,1);
    u_lid=X(i,2);
%      DisplayUV(RecU(:,:,i),RecV(:,:,i),RecP(:,:,i),Re,tf,dt,lx,ly,nx,ny,u_lid )

    TmepU=zeros(nxU+2,nyU+2);
    TmepU(:,end)=ones(nxU+2,1)*u_lid;
    TmepU(2:end-1,2:end-1)=RecU(:,:,i);
    U(i,:)=(TmepU(:))';
    
    TmepV=zeros(nxV+2,nyV+2);
%     TmepV(:,end)=ones(nV+2,1);
    TmepV(2:end-1,2:end-1)=RecV(:,:,i);
    V(i,:)=(TmepV(:))';

    TmepP=RecP(:,:,i);
    P(i,:)=(TmepP(:))';
     
end

options.nU=nxU+2;
options.dU=nyU+2;
options.nV=nxV+2;
options.dV=nyV+2;
options.nP=nxP;
options.dP=nyP;




end

