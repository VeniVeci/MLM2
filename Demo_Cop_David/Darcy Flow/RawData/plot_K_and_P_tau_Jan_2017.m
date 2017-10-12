
% Last modified: 25 Jan 2017 


% -----------------------   Plot permeability  and pressure fields ----------------------

perme_num = 20;  % Permeability and pressure fields to plot

% True fields
% load('./ZZZ_Third_Paper/data_lambda_0_3_sigma_1_0/K_true.mat'); K=Ktrue;  % Load permeabilities
% load('./ZZZ_Third_Paper/data_lambda_0_3_sigma_1_0/P_true.mat'); P=Ptrue;  % Load pressures

% Design fields
load('./RawData/data_lambda_0_3_sigma_1_0/K_design.mat'); K=Kdesign;  % Load permeabilities
load('./RawData/data_lambda_0_3_sigma_1_0/P_design.mat'); P=Pdesign;  % Load pressures

Nx=50;Ny=50;  % nodes each direction
Dx=1;Dy=1;      % [0,1]x[0,1]

% Plot Permeability field
Kmain=reshape(K(:,perme_num),50,50);

figure(1);
pcolor([0:Dx/(Nx-1):Dx],[0:Dy/(Ny-1):Dy],Kmain);
shading interp;
%shading flat;
%caxis([0.0 10.0]); % min(min(perme_field))  and max(max(perme_field))
colormap jet;
colorbar('FontSize',14);
xlabel('$x$','FontSize',18,'interpreter','latex');
ylabel('$z$','FontSize',18,'interpreter','latex');


% Plot corresponding pressure field for the given (above) permeability
Pmain=reshape(P(:,perme_num),50,50);

figure(2);
pcolor([0:Dx/(Nx-1):Dx],[0:Dy/(Ny-1):Dy],Pmain');
shading interp;
%shading flat;
%caxis([0.0 10.0]); % min(min(perme_field))  and max(max(perme_field))
colormap jet;
colorbar('FontSize',14);
xlabel('$x$','FontSize',18,'interpreter','latex');
ylabel('$z$','FontSize',18,'interpreter','latex');


% *************************** END CODE ***********************************
























