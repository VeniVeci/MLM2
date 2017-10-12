function [dz]=Burgers1D_DBC_FEM_ODE_MOR_DEIM_func(t,z,Mr,Br,Cr,Fr,v,U,UF,P)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %F = RHS F(f, in t);
    % For the moment assume F=0;

%     dy = M\ (F-(1/2)*B*y.^2-v*C*y);

%     Mr=U'*M*U;
%     Fr=U'*F;
%     Cr=U'*C*U;
%     Br=U*B;

    D=UF*inv(P'*UF);
    F_term_approx=D*(P'*(U*z).^2);
     
    dz = Mr\ (Fr-(1/2)*Br*F_term_approx-v*Cr*z);

%     dz = Mr\ (Fr-(1/2)*Br*(U*z).^2-v*Cr*z);
    
    

    
end