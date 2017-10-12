function [ Theata ] = Solver2D_MKII()
% 2D Richards equation with constant boundary condition 
% theata based Richards equation
%
% First edition: Weix 11/04/2017 
%%
tic
% nNode=10;
% nTime=100;

% Spatial setup
lengthZ=40;
deltaZ=5;
nNodeZ=lengthZ/deltaZ+1;

lengthX=40;
deltaX=5;
nNodeX=lengthX/deltaZ+1;

% Mesh
[x,z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);
indexMatrix = reshape(uint32(1:nNodeZ*nNodeX), nNodeZ, nNodeX);

% Temporal setup
lengthTime=1;
deltaTime=0.01;
nTime=lengthTime/deltaTime+1;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.01;



Theata=zeros(nNodeZ,nNodeX,nTime);
% H_temp=zeros(nNode,nMaxIteration);

%%% 
%initial state

%     Theata_init=zeros(nNodeZ,nNodeX);
%     Theata_init(:)=-61.5;
%     bcLeft=nNodeZ*ones(nNodeZ,1)*20.7;
%     bcRight=nNodeZ*ones(nNodeZ,1)*-61.5;
%     bcTop=nNodeX*ones(nNodeX,1)*20.7;
%     bcBottom=nNodeX*ones(nNodeX,1)*-61.5;
% 
%     % update initial field
%     Theata_init(1,:)=bcTop;
%     Theata_init(end,:)=bcBottom;
%     Theata_init(:,1)=bcLeft;
%     Theata_init(:,end)=bcRight;

theataS=0.48;
theata0=0.16;
mX=1;
mZ=1;

Theata_init=zeros(nNodeZ,nNodeX);

for j=1:nNodeX
    for k=1:nNodeZ
       xCoordi=deltaX*(j-1);
       zCoordi=deltaZ*(k-1);
       
       if (xCoordi/mX+zCoordi/mZ)<=1
           Theata_init(k,j)=  (theata0-theataS)*xCoordi/mX ...
                            + (theata0-theataS)*zCoordi/mZ ...
                            + theata0;
       else
           Theata_init(k,j)=theataS;
       end

    end
end


%% MAIN
TheataRecord(:,:,1)=Theata_init;
Theata=Theata_init;

for t=1:nTime
    
    TheataPreviousTime=Theata;
    
    K=K_Func(TheataPreviousTime);
    D=D_Func(TheataPreviousTime);

    index=0;

    for j=1:nNodeX
        for k=1:nNodeZ
            
            index=index+1; %next grid point 
            
%             %Coordinate to index
%             indexNextX=index+1;
%             indexLastX=index-1;
%             
%             indexNextZ=index+nNodeX+1;
%             indexLastZ=index-nNodeX-1;
            
            
            %--
            %check if top BC
            if k==1
                D_HalfLastZ=0;
                K_HalfLastZ=0;
                indexLastZ =0;
            else 
                D_HalfLastZ=(D(k,j)+D(k-1,j))/2; 
                K_HalfLastZ=(K(k,j)+K(k-1,j))/2;
                indexLastZ =index-1;
            end

            %check if bottom BC
            if k==nNodeZ
                D_HalfNextZ=0;
                K_HalfNextZ=0;
                indexNextZ =0;
            else 
                D_HalfNextZ=(D(k,j)+D(k+1,j))/2; 
                K_HalfNextZ=(K(k,j)+K(k+1,j))/2;
                indexNextZ =index+1;
            end

            %check if left BC
            if j==1
                D_HalfLastX=0;
                indexLastX =0;
            else 
                D_HalfLastX=(D(k,j)+D(k,j-1))/2; 
                indexLastX =index-nNodeZ;
            end

            %check if right BC
            if j==nNodeX
                D_HalfNextX=0;
                indexNextX=0;
            else 
                D_HalfNextX=(D(k,j)+D(k,j+1))/2; 
                indexNextX =index+nNodeZ;
            end       
            
            
%                 %--
%                 %check if top BC
%                 if indexLastZ<=0
%                     D_HalfLastZ=0;
%                     K_HalfLastZ=0;
%                 else 
%                     D_HalfLastZ=(D(k,j)+D(k-1,j))/2; 
%                     K_HalfLastZ=(K(k,j)+K(k-1,j))/2;
%                 end
% 
%                 %check if bottom BC
%                 if indexNextZ>=nNodeX*nNodeZ
%                     D_HalfNextZ=0;
%                     K_HalfNextZ=0;
%                 else 
%                     D_HalfNextZ=(D(k,j)+D(k+1,j))/2; 
%                     K_HalfNextZ=(K(k,j)+K(k+1,j))/2;
%                 end
% 
%                 %check if left BC
%                 if indexLastX<=0
%                     D_HalfLastX=0;
%                 else 
%                     D_HalfLastX=(D(k,j)+D(k-1,j))/2; 
%                 end
% 
%                 %check if right BC
%                 if indexNextX>=nNodeX*nNodeZ
%                     D_HalfNextX=0;
%                 else 
%                     D_HalfNextX=(D(k,j)+D(k+1,j))/2; 
%                 end                    

            %Assemble system 
            weightNextX=-deltaTime/deltaX^2*D_HalfNextX;
            weightLastX=-deltaTime/deltaX^2*D_HalfLastX;
            
            weightNestZ=-deltaTime/deltaZ^2*D_HalfNextZ;
            weightLastZ=-deltaTime/deltaZ^2*D_HalfLastZ;
            
            weightSelf = 1 - weightNextX ...
                           - weightLastX ...
                           - weightNestZ ...
                           - weightLastZ ;
            
                      
            b(index)=TheataPreviousTime(k,j)-deltaTime/deltaZ*(K_HalfNextZ-K_HalfLastZ);
            
            if indexNextX>0 A(index,indexNextX)=weightNextX; end
            if indexLastX>0 A(index,indexLastX)=weightLastX; end
            if indexNextZ>0 A(index,indexNextZ)=weightNestZ; end
            if indexLastZ>0 A(index,indexLastZ)=weightLastZ; end

        end

    end
    
    theata=A/b;
    Theata=reshape(theata,nNodeZ,nNodeX);
    
    TheataRecord(:,:,t)=Theata;
    
end
    
toc
    
figure
surf(Theata_init);

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



function result=theata(h)
theataS=0.287;
theataR=0.075;
alpha=1.611e6;
beta=3.96;

result=alpha*(theataS-theataR)/(alpha+abs(h)^beta)+theataR;
end

function result=theataDif(h)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

result=-alpha*(theata_s-theata_r)*-1*(alpha+abs(h)^beta)^(-2)*abs(h)^(beta-1);

end

function result=k(h)
rho=1.175e6;
r=4.74;
k_s=0.00944;

result=k_s*rho/(rho+abs(h)^r);
end