% Rec data process
%


clear 

load('exp2')

[nU,dU]=size(RecU(:,:,1));
[nV,dV]=size(RecV(:,:,1));

clearvars U V P

for i = 1:500
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


clearvars -except X U V P
% save('NS_Lid_Experiment1')