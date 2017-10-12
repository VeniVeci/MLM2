%% uxy_Demo

Paras.ampx=10;
Paras.ampy=10;

Paras.x0=0.5;
Paras.y0=0.5;


[X,Y] = meshgrid(0:0.1:1);

[Ux,Uy] = uxy(X,Y,Paras);
figure
quiver(X,Y,Ux,Uy)