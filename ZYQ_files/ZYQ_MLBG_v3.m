clear;
grPeriod_lpermm=2400;
grazing_angle_deg=2;

grBA_deg=1;
grAntiBA_deg=3;
photonEnergy_eV=2500;

material_sub='Si';
material_HZ='Cr';
material_LZ='C';

ML_d_nm=6.5;
ML_d_HZtod=0.45;
ML_N=60;

z_resolution_nm=0.1;
x_resolution_nm=1;
FourierOrders=5; %or FourierOrders=GR_order*GR_groove+3;
plSelect=0;

%estimate Inc
ML_order=1;
GR_groove=1;
GR_order=1;
thetaEst=estimateTheta(material_HZ,material_LZ, ...
    ML_order, ML_d_nm,ML_d_HZtod, ...
    GR_order,GR_groove,grPeriod_lpermm,photonEnergy_eV);

eff=[];
energy_=[];
effAll=[];
i=1;
outputData=zeros(length((thetaEst-0.2:0.01:thetaEst+0.2)),FourierOrders*2+2+1);

for grazing_angle_deg=(thetaEst-0.2:0.01:thetaEst+0.2)
ef = efficiency_bgrML(grPeriod_lpermm, ...
                grBA_deg, grAntiBA_deg, photonEnergy_eV, grazing_angle_deg, ...
                material_sub,material_HZ,material_LZ, ...
                ML_d_nm,ML_d_HZtod,ML_N,...
                z_resolution_nm,x_resolution_nm, FourierOrders,plSelect)
            outputData(i,1)=grazing_angle_deg;
            outputData(i,((1:FourierOrders)+1)*2)=ef.inc_top_reflected.efficiency;
            outputData(i,((1:FourierOrders)+1)*2-1)=ef.inc_top_reflected.theta;

eff(i,:)=eta;
theta(i)=grazing_angle_deg;
plot(theta,eff(:,1),theta,eff(:,2),'--',theta,eff(:,3),'*');xlabel('theta');title('Diffraction efficiency');legend('TE','TM','AVG');ylabel('-1th diffraction efficiency');pause(eps);

i=i+1;
end
% for grPeriod_lpermm=150 %[1000,500,150:-20:20]
%     effAll=[];
%     figure
%     for energy=500
%         for FourierOrders=1:50
%             photonEnergy_eV= energy;
%             eta = grating.efficiency_bgrML(grPeriod_lpermm, ...
%                 grBA_deg, grAntiBA_deg, photonEnergy_eV, grazing_angle_deg, ...
%                 material_sub,material_HZ,material_LZ, ...
%                 ML_d_nm,ML_d_HZtod,ML_N,...
%                 z_resolution_nm,x_resolution_nm, FourierOrders,plSelect);
%             eff=[eff,eta];
%             %     energy_=[energy_,energy];
%             
%         end
%         plot(1:50,eff,'-','DisplayName',['number of energy: ',num2str(energy)]);xlabel('Harmonics');ylabel('1 order efficiency''ordre 0');grid;title('TEST RETICOLO');pause(eps);
%         hold on
%         effAll=[effAll;eff];
%         eff=[];
%         energy_=[];
%     end
%     legend show
%     effAll=[1:50;effAll];
% %     save(['bgrCmp_effAll_v3_',num2str(grPeriod_lpermm),'LinePerMm.mat'],'effAll');
% end
% % 
function thetaEst=estimateTheta(material_HZ,material_LZ, ...
    ML_order, ML_d_nm,ML_d_HZtod, ...
    GR_order,GR_groove, grPeriod_lpermm,photonEnergy_eV)
    
for i=2:3
    if i==2
        nfile=['n_',material_HZ,'_cxro.txt'];
    elseif i==3
        nfile=['n_',material_LZ,'_cxro.txt'];
    end
    if exist(nfile, 'file') == 2
        nData=importdata(nfile);
        nData=nData.data;
    else
        disp('index File does not exist./type wrong material');
    end
    
    %wavelength
    lambda_nm = 1239.8/photonEnergy_eV;
    
    %index
    
    if i==2
        n_HZ_real=interp1(nData(:,1),nData(:,2),photonEnergy_eV);
        n_HZ_imag=interp1(nData(:,1),nData(:,3),photonEnergy_eV);
        n_HZ = 1-n_HZ_real+n_HZ_imag*1i;
    elseif i==3
        n_LZ_real=interp1(nData(:,1),nData(:,2),photonEnergy_eV);
        n_LZ_imag=interp1(nData(:,1),nData(:,3),photonEnergy_eV);
        n_LZ = 1-n_LZ_real+n_LZ_imag*1i;
    end
    
end

nAvg=n_LZ_real*(1-ML_d_HZtod)+n_HZ_real*ML_d_HZtod;

ThetaBRAGG = asind(sqrt((ML_order*1.2398/(2*ML_d_nm*photonEnergy_eV/1000))^2+2*nAvg)) ;
LD=1000/grPeriod_lpermm*1000;
thetaEst=(ThetaBRAGG-asind(GR_order*GR_groove*1.2398/LD/(photonEnergy_eV/1000)/2/sind(ThetaBRAGG))) ;


end

function ef = efficiency_bgrML(grPeriod_lpermm, ...
                grBA_deg, grAntiBA_deg, photonEnergy_eV, grazing_angle_deg, ...
                material_sub,material_HZ,material_LZ, ...
                ML_d_nm,ML_d_HZtod,ML_N,...
                z_resolution_nm,x_resolution_nm, FourierOrders,plSelect)

addpath(genpath('RETICOLO V7'))
retio

%number of Fourier orders
nn = FourierOrders;
%material  sub,HZ,LZ
for i=1:3
    if i==1
        nfile=['n_',material_sub,'_cxro.txt'];
    elseif i==2
        nfile=['n_',material_HZ,'_cxro.txt'];
    elseif i==3
        nfile=['n_',material_LZ,'_cxro.txt'];
    end
    if exist(nfile, 'file') == 2
        nData=importdata(nfile);
        nData=nData.data;
    else
        disp('index File does not exist./type wrong material');
    end
    
    %wavelength
    lambda_nm = 1239.8/photonEnergy_eV;
    
    %index
    if i==1
        n_sub_real=interp1(nData(:,1),nData(:,2),photonEnergy_eV);
        n_sub_imag=interp1(nData(:,1),nData(:,3),photonEnergy_eV);
        n_sub = 1-n_sub_real+n_sub_imag*1i;
    elseif i==2
        n_HZ_real=interp1(nData(:,1),nData(:,2),photonEnergy_eV);
        n_HZ_imag=interp1(nData(:,1),nData(:,3),photonEnergy_eV);
        n_HZ = 1-n_HZ_real+n_HZ_imag*1i;
    elseif i==3
        n_LZ_real=interp1(nData(:,1),nData(:,2),photonEnergy_eV);
        n_LZ_imag=interp1(nData(:,1),nData(:,3),photonEnergy_eV);
        n_LZ = 1-n_LZ_real+n_LZ_imag*1i;
    end
    
end
n_inc = 1;

%grating pitch
p_nm = 1/grPeriod_lpermm*1E6;

%angle of incidence
theta0_deg = 90-grazing_angle_deg;
k_parallel = n_inc*sin(theta0_deg*pi/180);

% thicknes mean grooveHeight
th_nm = p_nm*tand(grBA_deg)*tand(grAntiBA_deg)/(tand(grBA_deg)+tand(grAntiBA_deg));
th_nm = th_nm+z_resolution_nm*25+ML_d_nm*ML_N;

% generate texture and Profile
Apex_z=p_nm*tand(grBA_deg)*tand(grAntiBA_deg)/(tand(grBA_deg)+tand(grAntiBA_deg));
Apex_x=Apex_z/tand(grBA_deg);
Prf=[0,0;Apex_x,Apex_z;p_nm,0];

x=linspace(0,p_nm,round(p_nm/x_resolution_nm)+1);%
z=linspace(th_nm,0,round(th_nm/z_resolution_nm)+1);%

% genergy grid 
[X,Z]=meshgrid(x,z);

% number of blazed layers
N_layers=length(z);

% mat to store n 
n=X.*0;

% interplot the z of profilr based on x
Prf_z=interp1(Prf(:,1)',Prf(:,2)',x);

for i=1:ML_N
    Prf0=Prf_z+1+ML_d_nm*(i-1);
    Prf1=Prf_z+1+ML_d_nm*(i-1)+ML_d_nm*ML_d_HZtod;
    Prf2=Prf_z+1+ML_d_nm*i;
    if i==1
        P=find(Z<Prf0);         n(P)=n_sub;%substrate
        P=find(Z>=Prf0&Z<Prf1);  n(P)=n_HZ; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  n(P)=n_LZ; % spacing layer
    elseif i==ML_N %top layer
        P=find(Z>=Prf2);         n(P)=n_inc;% background
        P=find(Z>=Prf0&Z<Prf1);  n(P)=n_HZ; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  n(P)=n_LZ; % spacing layer     
    else
        P=find(Z>=Prf0&Z<Prf1);  n(P)=n_HZ; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  n(P)=n_LZ; % spacing layer     
    end
end

deltan = diff(n,1,2);
[nonZeroRowIndices, nonZeroColIndices]=find(deltan~= 0);
ProfEdge=[nonZeroRowIndices,nonZeroColIndices];
ProfEdge=sortrows(ProfEdge, 1);

textures = cell(1,N_layers+2);
textures{1}={n_inc};
textures{end}={n_sub};

for layer=1:N_layers
    if ~isnan(find(ProfEdge(:,1)==layer))
        p=find(ProfEdge(:,1)==layer);
        x_position=X(layer,ProfEdge(p,2)+1);
        n_value=n(layer,ProfEdge(p,2));     
        textures{layer+1} = {x_position, n_value};
    else
        n_value=n(layer,1);
        textures{layer+1} = {n_value};
    end
end


texture_list=1:N_layers+2;
th_list=[0,ones(1,N_layers).*z_resolution_nm,0];
profile={th_list,texture_list};% starting1, from incidence medium

% polarisation
i=1;
% for pol = [-1,1]; % 1:TE   -1:TM
pol = 1 ; % 1:TE   -1:TM
parm = res0(pol);  % parameter init

% initialisation
parm.res1.trace= 0;
aa = res1(lambda_nm, p_nm, textures, nn, k_parallel, parm);

ef = res2(aa, profile, parm);
% grating efficiency
% eta(i) = ef.inc_top_reflected.efficiency{-GR_groove*GR_order};

% i=i+1;
% end
% eta(3)=(eta(1)+eta(2))/2;

if plSelect==1
    %             show if needed
    % % xx=linspace(-p_nm,p_nm,length(x)*2);% on trace 2 periodes
    parm.res3.trace=1 ; % trace automatique
    parm.res3.cale=[];
    parm.res3.npts=[10,80,10];
    
    if pol==1 % 1:TE   -1:TM
        einc= ef.inc_top.PlaneWave_E(2);
    else
        einc= ef.inc_top.PlaneWave_H(2);
    end
    % x x section; aa from res1, profile, 1,
    [e,z,o]=res3(x,aa,profile,einc,parm);
    axis square
    set(gcf,'WindowStyle','docked')
else
end
% 
% 
 end