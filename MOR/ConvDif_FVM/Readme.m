%% Readme

%Simulating the unsteady state heat flow Diffusion Convection with velocity
%field specified 2D equation  
...by the Finite Volume Method using centreal differencing scheme in space
...and fully implicit scheme in time. Using the offered parameter.
    
% Governning equation : d(rho*T)/dt+k*Txx+k*Tyy=(rho*u*phi)x+(rho*u*phi)y+S; 
...Txx=dT/d^2(x);Tyy=dT/d^2(y);d(rho*u*phi)/dx;d(rho*u*phi)/dy;
...u=(u,y)' a vector indicates velocity.

% The point are discritized as follow grid. e.g a 3*4 grid. Numbers are
% index of point
% 
%        4--8--12
%        |  |   | 
%        3--7--11
%        |  |  |
%     ^  2--6--10
%     |  |  |  |
%  u_y|  1--5--9
%        ----> u_x
%
% Pe_x=F_x/D_x=(rho*u_x*dy)/(k*dy/dx)=rho*u_x*dx/k 
... Only if Pe_x<2 the result is stable and accurate in x direction
% Pe_y=F_y/D_y=(rho*u_y*dx)/(k*dx/dy)=rho*u_y*dy/k
... Only if Pe_y<2 the result is stable and accurate in y direction
%     
% Modifications:
% 6-April-2016, WeiX, first edition 