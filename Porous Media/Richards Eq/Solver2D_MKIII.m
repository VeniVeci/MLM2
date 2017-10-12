function [ H ] = Solver2D_MKIII()
% 2D Richards equation with constant boundary condition 
% h based Richards equation
%
% First edition: Weix 11/04/2017 
%%
tic
% nNode=10;
% nTime=100;

% Spatial setup
lengthZ=40;
deltaZ=8;
nNodeZ=lengthZ/deltaZ+1;

lengthX=40;
deltaX=5;
nNodeX=lengthX/deltaX+1;

% Mesh
[x,z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);
indexMatrix = reshape(uint32(1:nNodeZ*nNodeX), nNodeZ, nNodeX);

% Temporal setup
lengthTime=10;
deltaTime=1;
nTime=lengthTime/deltaTime+1;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.01;


H=zeros(nNodeZ,nNodeX,nTime);
% H_temp=zeros(nNode,nMaxIteration);

%%% 
%initial state

H_init=zeros(nNodeZ,nNodeX);
H_init(:)=-61.5;
bcLeft=nNodeZ*ones(nNodeZ,1)*20.7;
bcRight=nNodeZ*ones(nNodeZ,1)*-61.5;
bcTop=nNodeX*ones(nNodeX,1)*20.7;
bcBottom=nNodeX*ones(nNodeX,1)*-61.5;

% update initial field
H_init(1,:)=bcTop;
H_init(end,:)=bcBottom;
H_init(:,1)=bcLeft;
H_init(:,end)=bcRight;






%% MAIN
TheataRecord(:,:,1)=H_init;
H=H_init;

% C=ones(nNodeZ,nNodeX)*1234567;

for t=1:nTime
    
    H_PreviousTime=H;
    
    
    for k=1:nMaxIteration 
    
        H0=H;
        
        C=theataDifFunc(H);
        K=kFunc(H);
        D=D_Func(H);


        %initialize
        A=zeros(nNodeZ*nNodeX);
        B=zeros(nNodeZ*nNodeX,1);

        %Assemble Ax+B=0;
        for i =2:size(indexMatrix,2)-1
            for j =2:size(indexMatrix,1)-1

                kHalfUp   =(K(j,i)+K(j-1,i))/2;
                kHalfDown =(K(j,i)+K(j+1,i))/2;
                kHalfLeft =(K(j,i)+K(j,i-1))/2;
                kHalfRight=(K(j,i)+K(j,i+1))/2;

                cCenter=C(j,i);


                wUp   = -kHalfUp  ./deltaZ^2;
                wDown = -kHalfDown./deltaZ^2;
                wLeft = -kHalfLeft./deltaX^2;
                wRight= -kHalfRight./deltaX^2;

                wCenter=cCenter/deltaTime-wUp-wDown-wLeft-wRight;

                b=(kHalfDown-kHalfUp)/deltaZ-H_PreviousTime(j,i)*cCenter/deltaTime;


                indexCenter=indexMatrix(j,i);
                indexUp=indexMatrix(j-1,i);
                indexDown=indexMatrix(j+1,i);
                indexLeft=indexMatrix(j,i-1);
                indexRight=indexMatrix(j,i+1);

                %Check BC and modify 
                %check if top BC
                if j==2 
                    indexUp=0;
                    b=b+wUp*H(j-1,i);
                end

                %check if bottom BC
                if j==nNodeZ-1 
                    indexDown=0;
                    b=b+wDown*H(j+1,i);
                end

                %check if left BC
                if i==2 
                    indexLeft=0;
                    b=b+wLeft*H(j,i-1);
                end

                %check if right BC
                if i==nNodeX-1 
                    indexRight=0;
                    b=b-wRight*H(j,i+1);
                end


                if indexUp>0 A(indexCenter,indexUp)=wUp; end
                if indexDown>0 A(indexCenter,indexDown)=wDown; end
                if indexLeft>0 A(indexCenter,indexLeft)=wLeft; end
                if indexRight>0 A(indexCenter,indexRight)=wRight; end

                B(indexCenter,1)=b;


            end
        end
        %trime for BC
        indexIfKnown=~any(A);
        indexIfUnknown=any(A);

        A=A(indexIfUnknown,indexIfUnknown);
        B=B(indexIfUnknown);


        h=A\-B;
        H(2:end-1,2:end-1)=reshape(h,nNodeZ-2,nNodeX-2);
    
        
        sseIte=sum((H(:)-H0(:)).^2);
        if sqrt(sseIte)<miniIteError 
            break 
        end
        
    end
    
    
    
    
    TheataRecord(:,:,t)=H;
    
end
    
toc
    
figure
surf(H_init);

for t=1:nTime
    surf(TheataRecord(:,:,t))
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
    
end
    

end


function [k,j]=index2kj(index)
    


end



function K=K_Func(Theata)
b=1;
K_S=00944;
% Phi_S=1;
theataS=0.48;

K=K_S.*(Theata./theataS).^(2.*b+3);

end

function D=D_Func(Theata)
b=1;
K_S=00944;
Phi_S=1;
theata_S=0.48;

D=-(b.*K_S.*Phi_S./theata_S)*(Theata./theata_S).^(b+2);

end



function theata=theataFunc(h)
theataS=0.287;
theataR=0.075;
alpha=1.611e6;
beta=3.96;

result=alpha.*(theataS-theataR)/(alpha+abs(h).^beta)+theataR;
end

function theataDif=theataDifFunc(h)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

theataDif=-alpha.*(theata_s-theata_r).*-1.*(alpha+abs(h).^beta).^(-2).*abs(h).^(beta-1);

end

function result=kFunc(h)
rho=1.175e6;
r=4.74;
k_s=0.00944;

result=k_s.*rho./(rho+abs(h).^r);
end