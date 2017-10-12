function [ H ] = Solver_MKI_fminunc()
% 1D Richards equation with constant boundary condition 
%
% First edition: Weix 23/03/2017 
%%
tic
nNode=10;
nTime=100;

% Spatial setup
lengthZ=40;
deltaZ= 8;
nNode=lengthZ/deltaZ+1;

% Temporal setup
lengthTime=360;
deltaTime=1;
nTime=lengthTime/deltaTime+1;

% Iteration setup
nMaxIteration=1000;
miniIteError=0.01;



H=zeros(nNode,nTime);
H_temp=zeros(nNode,nMaxIteration);

%%% 

% H=zeros(nTime,nNode,nMaxIteration);
% time=[0:deltaTime:100]

% Initial condition and Boundary condition
h_init=zeros(nNode,1);
h_init(:)=-61.5;
h_init(1)=-20.7;
h_init(end)=-61.5;


%% MAIN

hTop=-20.7;
hEnd=-61.5;
h=zeros(nNode-2,1);
h(:)=-61.5;

for t=1:nTime
    
    h0=h;
    fun = @(h) residualFunc(h,h0,hTop,hEnd,deltaZ,deltaTime);

%     options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
    
%     h = fzero(fun,h0); %for scalar value 
%     h = fsolve(fun,h0,options);

%     options = optimset('LargeScale','off', 'HessUpdate','bfgs' );
%     options = optimset(options,'gradobj','on');
    
    options.OptimalityTolerance=1e-3;
    h = fminunc(fun,h0,options);
     
    H(:,t+1)=[hTop;h;hEnd];
    
end
computationalTime=toc
    
figure

for t=1:nTime
    plot(H(:,t))
    title(sprintf('time=%i',t))
    drawnow
    frame(t)=getframe;
    
end
    
   
end

function residual= residualFunc(h,h0,hTop,hEnd,deltaZ,deltaTime)
    r=theata(h)- deltaTime.*rightDif(h,hTop,hEnd,deltaZ)-deltaTime.*theata(h0);
    residual=r'*r;
    
end

function dif= rightDif(h,hTop,hEnd,deltaZ)

    h_all=[hTop;h;hEnd];
    
    for i =2:length(h_all)-1

        dif(i-1,1)=((k(h_all(i+1))+k(h_all(i)))/2*(h_all(i+1)-h_all(i))...
                -(k(h_all(i))+k(h_all(i-1)))/2*(h_all(i)-h_all(i-1)))/deltaZ^2 ...
                -(k(h_all(i+1))-k(h_all(i-1)))/(2*deltaZ);
    end
end



function result=theata(h)
theata_s=0.287;
theata_r=0.075;
alpha=1.611e6;
beta=3.96;

result=alpha*(theata_s-theata_r)./(alpha+abs(h).^beta)+theata_r;
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

result=k_s.*rho./(rho+abs(h).^r);
end