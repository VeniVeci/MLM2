
%------------------------------------------------------------------------
%--- Heat Equation in two dimensions-------------------------------------
%--- Solves Ut=alpha*(Uxx+Uyy)-------------------------------------------
%------------------------------------------------------------------------

clc;
close all;
clear all;

%--dimensions...........................................................

N = 40;  
DX=0.1; % step size
DY=0.1;
Nx=5; 
Ny=5;

X=0:DX:Nx;
Y=0:DY:Ny;



alpha=5; % arbitrary thermal diffusivity 

%--boundary conditions----------------------------------------------------

U(1:N,1:N) = 0 ;

U(1,1:N) = 100; 
U(N,1:N) = 300;  
U(1:N,1) = 200;  
U(1:N,N) = 400;  

%--initial condition------------------------------------------------------

U(23:29,23:29)=0; % a heated patch at the center

Umax=max(max(U));

%-------------------------------------------------------------------------

DT = DX^2/(2*alpha); % time step 

M=2000; % maximum number of allowed iteration

%---finite difference scheme----------------------------------------------

fram=0;
Ncount=0;
loop=1;
while loop==1;
   ERR=0; 
   U_old = U;
for i = 2:N-1
for j = 2:N-1
   Residue=(DT*((U_old(i+1,j)-2*U_old(i,j)+U_old(i-1,j))/DX^2 ... 
                      + (U_old(i,j+1)-2*U_old(i,j)+U_old(i,j-1))/DY^2) ...
                      + U_old(i,j))-U(i,j);
   ERR=ERR+abs(Residue);
  U(i,j)=U(i,j)+Residue;
  
  % U(i,j)=(DT*((U_old(i+1,j)-2*U_old(i,j)+U_old(i-1,j))/DX^2 + (U_old(i,j+1)-2*U_old(i,j)+U_old(i,j-1))/DY^2) + U_old(i,j))
 
end
end

if(ERR>=0.01*Umax)  % allowed error limit is 1% of maximum temperature
    Ncount=Ncount+1;
         if (mod(Ncount,50)==0) % displays movie frame every 50 time steps
              fram=fram+1;
              surf(U);
              axis([1 N 1 N ])
              h=gca; 
              get(h,'FontSize') 
              set(h,'FontSize',12)
              colorbar('location','eastoutside','fontsize',12);
              xlabel('X','fontSize',12);
              ylabel('Y','fontSize',12);
              title('Heat Diffusion','fontsize',12);
              fh = figure(1);
             set(fh, 'color', 'white'); 
            F=getframe;
         end
 
 %--if solution do not converge in 2000 time steps------------------------
 
    if(Ncount>M)
        loop=0;
        disp(['solution do not reach steady state in ',num2str(M),...
            'time steps'])
    end
    
 %--if solution converges within 2000 time steps..........................   
    
else
    loop=0;
    disp(['solution reaches steady state in ',num2str(Ncount) ,'time steps'])
end
end
%------------------------------------------------------------------------
%--display a movie of heat diffusion------------------------------------

 movie(F,fram,1)

%------END---------------------------------------------------------------




