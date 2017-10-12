function [dy]=Burgers1D_DBC_FEM_ODE_DEIM_func(t,y,M,B,C,F,v,Dr,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Burgers' equation 1D, finite element method, Dirichlet boundary 
% conditions, solver.
% 
% Function inputs:
%
% t : current time
% y : current y
% v : viscosity
%
%
% Function outputs:
% out y : computed derivative
%11-Sep-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %F = RHS F(f, in t);
    % For the moment assume F=0;
    
%     dy = M\ (F-(1/2)*B*y.^2-v*C*y);
%     
%     dy = M\ (F-(1/2)*(Dr*P'*(B*y.^2))-v*C*y);
    
    dy = M\ (F-(1/2)*B* (Dr*P'*(y.^2)) -v*C*y);
 
end