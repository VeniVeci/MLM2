function [] = deimTest()
%% 
% Test DEIM-POD method when variating the number of basis of DEIM and 
% POD mode
%
% 08/02/2017 Weix First edition 
%
%
%% Setup
% ------------------Problem Parameters------------------------------------- 
Paras.Re=1000;    % Reynolds Number
% v=1/Re;     % viscosity

% ------------------Solver Parameters--------------------------------------
Paras.n=128;           % Total Spatial elements
Paras.t_end=4;        % End time
Paras.t_n=100;        % Number of time step 
% Paras.t=0:(t_end/t_n):t_end; % time sequence (t=0 does not account 1 step)

% solver = 'ode45';
% options = odeset('RelTol',1e-6,'AbsTol',1e-10);

nPod=20:1:100;
nDeim=20:5:50;


%Initial condition
u0a=5; 
u0b=1; 
u0=@(x) 1*exp(-(u0a*x+u0b)).*sin(3*pi*x);
 
%Sources term
gx=@(x) 0.02*exp(x);


%% Main 

% Normal solver
[Y_Ture,T,TimeCost_FOM]=burgerSolver(Paras,gx,u0);

% POD on snapshots (use all available)
Y_inter=Y_Ture(2:end-1,:);
[U,S,~]=svd(Y_inter,'econ');  % U*S*V'=Rec_X
% U=U(:,1:approximate_degree);
eigenvalues=diag(S);


% POD on non-linear terms. This is only valid for the Burgers equation
% where the nonlinear term is element-wise square of u(x)
[U_DEIM,S_DEIM,~]=svd(Y_inter.^2); 
[~,U_DEIM,P] = DEIM(U_DEIM);


% running with different setting
nIteration=length(nPod)*length(nDeim);
h = waitbar(0,'');

for i=1:length(nPod)
    for j=1:length(nDeim)
                
        [Y_Predict,T3,TimeCost_DEIM_POD]=burgerPod_deimSolver(Paras,gx,u0,U(:,1:nPod(i)),U_DEIM(:,1:nDeim(j)),P(:,1:nDeim(j)));
        RE(i,j)=SpaceTimeRelativeError(Y_Predict,Y_Ture);
        
        % wait bar display
        iIteration=(i-1)*length(nDeim)+j;      
        string=sprintf('%i/%i iteration done',iIteration,nIteration);    
        waitbar(iIteration/nIteration,h,string);
        
    end
end
close(h)


nTick=10;

figure(1)
bar3(RE)
ylabel('nPod')
xlabel('nDeim')
zlabel('relative error')

ax = gca; 
% %   ax.ZScale='log';

ax.YTick=1:10:length(nPod);
ax.YTickLabel=nPod(1:10:end);

ax.XTickLabel=nDeim;




figure(2)
mesh(RE)
ylabel('nPod')
xlabel('nDeim')
zlabel('relative error')

ax = gca; 
ax.ZScale='log';
ax.YTickLabel=nPod;
ax.XTickLabel=nDeim;
end


function [Error]=SpaceTimeRelativeError(y,yTure)
% Y= dimensionality x time matrix
squareError=(y-yTure).^2;
squareSumError=sum(squareError);
relativeError=sqrt(squareSumError./sum(yTure.^2));
Error=mean(relativeError);

end 







