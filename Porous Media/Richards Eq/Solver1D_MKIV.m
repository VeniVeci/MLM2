function [  ] = Solver1D_MKIV()
% 1D Richards equation with constant boundary condition 
% h based Richards equation
%
% First edition: Weix 12/04/2017 
%%
tic
% nNode=10;
% nTime=100;

% Spatial setup
lengthZ=40;
deltaZ=1;
nNodeZ=lengthZ/deltaZ+1;

% Temporal setup
lengthTime=360;
deltaTime=1;
nTime=lengthTime/deltaTime+1;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.01;


% Mesh
% [x,z] = meshgrid(0:deltaX:lengthX,0:deltaZ:lengthZ);

% integer shows unknowns and 0 shows DBC
nUnknown=(nNodeZ-2);
indexUnknown=(1:nNodeZ-2);
% indexUnknown = reshape(uint32(1:nUnknown), (nNodeZ-2));

indexMatrix=zeros(nNodeZ,1);
indexMatrix(2:end-1)=indexUnknown;



H=zeros(nNodeZ,nTime);
% H_temp=zeros(nNode,nMaxIteration);

%%% 
%initial state

H_init=zeros(nNodeZ,1);
H_init(:)=-61.5;

bcTop=-20.7;
bcBottom=-61.5;
% bcBottom=-20.7;

% update initial field
H_init(1,:)=bcTop;
H_init(end,:)=bcBottom;





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


        %initialize
        A=zeros(nUnknown);
        B=zeros(nUnknown,1);

        %Assemble Ax+B=0;
            
            i=1;
            for j =2:size(indexMatrix,1)-1

                kHalfUp   =(K(j,i)+K(j-1,i))/2;
                kHalfDown =(K(j,i)+K(j+1,i))/2;

                cCenter=C(j,i);


                wUp   = -kHalfUp  ./deltaZ^2;
                wDown = -kHalfDown./deltaZ^2;

                wCenter=cCenter/deltaTime-wUp-wDown;

                b=(kHalfDown-kHalfUp)/deltaZ-H_PreviousTime(j,i)*cCenter/deltaTime;


                indexCenter=indexMatrix(j,i);
                indexUp=indexMatrix(j-1,i);
                indexDown=indexMatrix(j+1,i);


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


                if indexUp>0 A(indexCenter,indexUp)=wUp; end
                if indexDown>0 A(indexCenter,indexDown)=wDown; end

                A(indexCenter,indexCenter)=wCenter;
                B(indexCenter,1)=b;

            end
 
        
%         %trime for BC
%         indexIfKnown=~any(A);
%         indexIfUnknown=any(A);
% 
%         A=A(indexIfUnknown,indexIfUnknown);
%         B=B(indexIfUnknown);


        h=A\(-B);
        H(2:end-1)=h;
    
        
        sseIte=sum((H(:)-H0(:)).^2);
        if sqrt(sseIte)<miniIteError 
            break 
        end
        
    end
    
    
    
    
    TheataRecord(:,:,t)=H;
    
end
    
toc
    
figure
plot(H_init);

for t=1:nTime
    plot(TheataRecord(:,:,t))
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
    
end
    

end





function [k,j]=index2kj(index)
    
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

