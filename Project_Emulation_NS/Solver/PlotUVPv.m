function [  ] = PlotUVPv(U_in,V_in,P_in,Re,u_lid,options)
%Plot field graph using U V P (vector form)
% U V P is the form of vector 
% PlotUVPv(U_starorig(1,:),V_starorig(1,:),P_starorig(1,:),options,Name_field)

U_field = reshape(U_in,[options.nU,options.dU]);
V_field = reshape(V_in,[options.nV,options.dV]);
P_field = reshape(P_in,[options.nP,options.dP]);

U_field=U_field(1:end-2,1:end-2);
V_field=V_field(1:end-2,1:end-2);

DisplayUV( U_field,V_field,P_field,Re,200,options.dt,options.lx,options.ly,options.nx,options.ny,u_lid )
% constant 200 is due to a program bug. it should be 'tf' which would
% cause the bug.


end

