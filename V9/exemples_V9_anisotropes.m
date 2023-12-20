% exemples_V9_anisotropes
%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 1D (TE or TM) %
%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength=8;
period=10;% same unit as wavelength
n_incident_medium=1;% refractive index of the top layer
n_transmitted_medium=1.5;% refractive index of the bottom layer
angle_theta0=-10;k_parallel=n_incident_medium*sin(angle_theta0*pi/180);
parm=res0(1);parm.not_io=1;% TE polarization. For TM : parm=res0(-1);parm.not_io=1;
parm.res1.champ=1;% the electromagnetic field is calculated accurately
nn=40;% Fourier harmonics run from [-40,40]
% textures for all layers including the top and bottom layers
texture=cell(1,3);
textures{1}= n_incident_medium; % uniform texture
textures{2}= n_transmitted_medium; % uniform texture
textures{3}={[-2.5,2.5],[n_incident_medium,n_transmitted_medium] };
[aa,neff]=res1(wavelength,period,textures,nn,k_parallel,parm);
Num_texture=3;Num_mode=1;res1(aa,neff, Num_texture, Num_mode); % To plot the profile of Bloch mode Num_mode of the texture Num_texture


profile={[4.1,5.2,4.1],[1,3,2]};
one_D_TE=res2(aa,profile)
eff=one_D_TE.inc_top_reflected.efficiency{-1}
J=one_D_TE.Jones.inc_top_reflected{-1};% Jones’coefficients
abs(J)^2 % first order efficiency for an illumination from the top layer
% field calculation
x=linspace(-period/2,period/2,51);% x coordinates(z-coordinates are determined by res3.m)
einc=1;
parm.res3.trace=1; % plotting automatically
parm.res3.npts=[50,50,50];
[e,z,index]=res3(x,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,1)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ey)');

% Loss calculation
n_anisotropic=9.999;parm.res1.change_index={[n_anisotropic,.1+5i,.5+1i,2+5i]};
textures{3}={[-2.5,2.5],[n_incident_medium,n_anisotropic] };
aa_loss=res1(wavelength,period,textures,nn,k_parallel,parm);
one_D_loss=res2(aa_loss,profile)
parm.res3.npts=[[5,10,5];[3,3,3]];
einc=one_D_loss.inc_top.PlaneWave_E(2);
[e,z,index,wZ,loss_per_layer,loss_of_Z,loss_of_Z_X,X,wX]=res3([-period/2,period/2],aa_loss,profile,einc,parm);
Energie_conservation=sum(one_D_loss.inc_top_reflected.efficiency)+sum(one_D_loss.inc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5* period)-1
% Loss with Poynting
parm.res3.trace=0; 
parm.res3.pertes_poynting=1;
[e,z,index,wZ,loss_per_layer]=res3([-period/2,period/2],aa_loss,profile,einc,parm);
Energie_conservation_Poynting=sum(one_D_loss.inc_top_reflected.efficiency)+sum(one_D_loss.inc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5* period)-1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMPLE EXAMPLE 1D CONICAL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength=8;
period=10;% same unit as wavelength
n_incident_medium=1;%refractive index of the top layer
n_transmitted_medium=1.5;% refractive index of the bottom layer
angle_theta0=10;k_parallel=n_incident_medium*sin(angle_theta0*pi/180);
angle_delta=-20;
parm=res0;parm.not_io=1; % default parameters for "parm"
parm.res1.champ=1; % the electromagnetic field is calculated accurately
nn=5; % Fourier harmonics run from [-5,5]
% textures for all layers including the top and bottom layers
texture=cell(1,3);
textures{1}= n_incident_medium; % uniform texture
textures{2}= n_transmitted_medium; % uniform texture
textures{3}={[-2.5,2.5],[n_incident_medium,n_transmitted_medium] };
[aa,neff]=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);
Num_texture=3;Num_mode=1;res1(aa,neff, Num_texture, Num_mode,linspace(-period/2,period/2,100)); % To plot the profile of Bloch mode Num_mode of the texture Num_texture

profile={[4.1,5.2,4.1],[1,3,2]};
conical=res2(aa,profile)
eff_TETM=conical.TEinc_top_reflected.efficiency{-1}
% -1 order efficiency (TE+TM) for a TE-illumination from the top layer
eff_TE=conical.TEinc_bottom_transmitted.efficiency_TE{-1}
% -1 order TE efficiency for a TE-illumination from the top layer
J=conical.Jones.inc_bottom_transmitted{-1};% Jones’matrix
abs(J).^2 % -1 order efficiencies for an illumination from the top layer
% field calculation
x=linspace(-period/2,period/2,51);% x coordinates(z-coordinates are determined by res3.m)
einc=[0,1]; % E-field components in the (u, v) basis (default is illumination from the top layer)
parm.res3.trace=1; % plotting automatically
parm.res3.npts=[50,50,50];
[e,z,index]=res3(x,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,3)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ez)');

% Loss calculation
n_anisotropic=rand;parm.res1.change_index={[n_anisotropic,.1+5i,.5+1i,2+5i]};
textures{3}={[-2.5,2.5],[n_incident_medium,n_anisotropic] };
aa_loss=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);
conical_loss=res2(aa_loss,profile)
parm.res3.npts=[[5,10,5];[4,10,4]];
einc= conical_loss.TEinc_top.PlaneWave_TE_Eu;
[e,z,index,wZ,loss_per_layer,loss_of_Z,loss_of_Z_X,X,wX]=res3([-period/2,period/2],aa_loss,profile,einc,parm);
Energie_conservation=sum(conical_loss.TEinc_top_reflected.efficiency)+sum(conical_loss.TEinc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5* period)-1
% Loss with Poynting
parm.res3.trace=0; 
parm.res3.pertes_poynting=1;
[e,z,index,wZ,loss_per_layer]=res3([-period/2,period/2],aa_loss,profile,einc,parm);
Energie_conservation_Poynting=sum(conical_loss.TEinc_top_reflected.efficiency)+sum(conical_loss.TEinc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5* period)-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% SIMPLE EXAMPLE 2D %
%%%%%%%%%%%%%%%%%%%%%
wavelength=8;
period=[10,15];% same unit as wavelength
n_incident_medium=1;% refractive index of the top layer
n_transmitted_medium=1.5;% refractive index of the bottom layer
angle_theta=10;k_parallel=n_incident_medium*sin(angle_theta*pi/180);
angle_delta=-20;
parm=res0;parm.not_io=1; % default parameters for "parm"
parm.res1.champ=1; % the eletromagnetic field is calculated accurately
nn=[3,2]; % Fourier harmonics run from [-3,3]in x and [-2,2] in y
% textures for all layers including the top and bottom
texture=cell(1,3);
textures{1}= n_incident_medium; % uniform texture
textures{2}= n_transmitted_medium; % uniform texture
textures{3}={n_incident_medium,[0,0,5,2,n_transmitted_medium,1] };
[aa,neff]=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);
Num_texture=3;Num_mode=1;res1(aa,neff, Num_texture, Num_mode); % To plot the profile of Bloch mode Num_mode of the texture Num_texture


profile={[4.1,5.2,4.1],[1,3,2]};
two_D=res2(aa,profile)
eff_TETM=two_D.TEinc_top_reflected.efficiency{-1,1}
% (-1,1) order efficiency (TE+TM) for a TE-illumination from the top layer
eff_TE=two_D.TEinc_bottom_transmitted.efficiency_TE{-1,1}
% (-1,1) TE efficiency for a TE-illumination from the top layer
J=two_D.Jones.inc_bottom_transmitted{-1,1};% Jones’matrix
abs(J).^2 % (-1,1) order efficiency for an illumination from the bottom layer
% field calculation in plane y=0
x=linspace(-period(1)/2,period(1)/2,51);y=0;%(x,y) coordinates (z-coordinates are determined by res3.m)
einc=[0,1];% E-field components in the (u, v) basis (default is illumination from the top layer)
parm.res3.trace=1; % plotting automatically
parm.res3.npts=[50,50,50];
[e,z,index]=res3(x,y,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,:,2)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ey)');
% Loss calculation
n_anisotropic=rand;parm.res1.change_index={[n_anisotropic,.1+5i,.5+1i,2+5i]};
textures{3}={n_anisotropic,[0,0,5,2,1,1] };
aa_loss=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);
two_D_loss=res2(aa_loss,profile)
parm.res3.npts=[[0,10,0];[1,3,1]];
einc= two_D_loss.TEinc_top.PlaneWave_TE_Eu;
[e,z,index,wZ,loss_per_layer,loss_of_Z,loss_of_Z_X_Y,X,Y,wXY]=res3([-period(1)/2,period(1)/2],[-period(2)/2,period(2)/2],aa_loss,profile,einc,parm);
Energie_conservation=sum(two_D_loss.TEinc_top_reflected.efficiency)+sum(two_D_loss.TEinc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5*prod(period))-1
% Loss with Poynting
parm.res3.trace=0; 
parm.res3.pertes_poynting=1;
[e,z,index,wZ,loss_per_layer]=res3([-period(1)/2,period(1)/2],[-period(2)/2,period(2)/2],aa_loss,profile,einc,parm);
Energie_conservation_Poynting=sum(two_D_loss.TEinc_top_reflected.efficiency)+sum(two_D_loss.TEinc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5*prod(period))-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMPLE EXAMPLE 0D CONICAL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength=8;
period=10;% same unit as wavelength
n_incident_medium=1;%refractive index of the top layer
n_transmitted_medium=1.5;% refractive index of the bottom layer
angle_theta0=10;k_parallel=n_incident_medium*sin(angle_theta0*pi/180);
angle_delta=-20;
parm=res0;parm.not_io=1; % default parameters for "parm"
parm.res1.champ=1; % the electromagnetic field is calculated accurately
nn=0; % Fourier harmonics only 0
% textures for all layers including the top and bottom layers
texture=cell(1,3);
textures{1}= n_incident_medium; % uniform texture
textures{2}= n_transmitted_medium; % uniform texture
epsilon=[[2.1160         0    0.7165];[0    1.3995         0]; [0.7165         0    2.1160]];
textures{3}={ epsilon} ;
[aa,neff]=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);




profile={[4.1,5.2,4.1],[1,3,2]};
conical=res2(aa,profile)
% field calculation
x=linspace(-period/2,period/2,51);% x coordinates(z-coordinates are determined by res3.m)
einc=[0,1];
 % E-field components in the (u, v) basis (default is illumination from the top layer)
parm.res3.trace=1; % plotting automatically
parm.res3.npts=[50,50,50];
[e,z,index]=res3(x,aa,profile,einc,parm);
figure;pcolor(x,z,real(squeeze(e(:,:,3)))); % user plotting
shading flat;xlabel('x');ylabel('y');axis equal;title('Real(Ez)');

% Loss calculation
epsilon=randn(3)+1i*randn(3);epsilon=epsilon+epsilon';H=randn(3,1)+1i*randn(3,1);epsilon=1i*H*H'+epsilon';
% anisotropie generale non amplificatrice
textures{3}={epsilon };
%textures{3}=1+5i ;
aa_loss=res1(wavelength,period,textures,nn,k_parallel,angle_delta,parm);
conical_loss=res2(aa_loss,profile)
einc=conical_loss.TEinc_top.PlaneWave_TE_Eu;
parm.res3.npts=[[5,10,5];[4,10,4]];

% Loss with Poynting
parm.res3.trace=0; 
parm.res3.pertes_poynting=1;
[e,z,index,wZ,loss_per_layer]=res3([-period/2,period/2],aa_loss,profile,einc,parm);
Energie_conservation_Poynting=sum(conical_loss.TEinc_top_reflected.efficiency)+sum(conical_loss.TEinc_top_transmitted.efficiency)+sum(loss_per_layer)/(.5* period)-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%