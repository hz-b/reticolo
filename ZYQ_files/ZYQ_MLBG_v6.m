clear;

addpath(genpath(fullfile(pwd, 'helpers')))

metadata.grPeriod_lpermm=2400;
metadata.grBA_deg=1;
metadata.grAntiBA_deg=3;

metadata.material_sub='Si';
metadata.material_HZ='Cr';
metadata.material_LZ='C';

metadata.ML_d_nm=6.5;
metadata.ML_d_HZtod=0.45;
metadata.ML_N=60;

metadata.z_resolution_nm=0.01;
metadata.x_resolution_nm=1;
metadata.FourierOrders=5; %or FourierOrders=GR_order*GR_groove+3;
metadata.plSelect=1;

%estimate Inc
metadata.ML_order=1;
metadata.GR_groove=1;
metadata.GR_order=1;
metadata.photonEnergy_eV=2500;


% Find parameters with more than one value
fields = fieldnames(metadata);
disp('Parameters with more than one value:');
loopNumber=1;
for i = 1:numel(fields)
    field = fields{i};
    if ~ischar(metadata.(field))&& numel(metadata.(field)) > 1
        loop.para{loopNumber}=field;
        loop.Idx{loopNumber}=i;
        loop.value{loopNumber}=metadata.(field);
        disp(['Field: ' field]);
        disp(['Values: ' num2str(metadata.(field))]);
        loopNumber=loopNumber+1;
    end
end
loopNumber=loopNumber-1;

% Define loop bounds based on the given number
loopBounds = zeros(1, loopNumber);
for i = 1:loopNumber
    IdxCurrent=loop.Idx{loopNumber};
    loopBounds(i) = length(metadata.(fields{IdxCurrent}));
end

%outputdata 0 to store all data, 1 to store max eff
outputdata0=cell(1,prod(loopBounds, 2));

% Generate nested loops based on the specified bounds
numLoops = length(loopBounds);
outputdata1=zeros(prod(loopBounds, 2),numLoops+3);

% Initialize cell array to store loop indices
loopIndices = cell(1, numLoops);

% Initialize loop indices
for i = 1:numLoops
    loopIndices{i} = 1;
end

% Nested loop structure
k=1;
idx_outputdata0=1;

% Start the timer
tic;
if numLoops==0
    disp('no Loop indices, single point calculation ');
    metadataCurrent=metadata;
       thetaEst=estimateTheta(metadataCurrent.material_HZ,metadataCurrent.material_LZ, ...
        metadataCurrent.ML_order, metadataCurrent.ML_d_nm,metadataCurrent.ML_d_HZtod, ...
        metadataCurrent.GR_order,metadataCurrent.GR_groove,metadataCurrent.grPeriod_lpermm,metadataCurrent.photonEnergy_eV);

    grazing_angle_deg=thetaEst; %(thetaEst-0.06*2:0.005:thetaEst+0.06*2)
        ef = efficiency_bgrML(metadataCurrent.grPeriod_lpermm, ...
            metadataCurrent.grBA_deg, metadataCurrent.grAntiBA_deg, metadataCurrent.photonEnergy_eV, grazing_angle_deg, ...
            metadataCurrent.material_sub,metadataCurrent.material_HZ,metadataCurrent.material_LZ, ...
            metadataCurrent.ML_d_nm,metadataCurrent.ML_d_HZtod,metadataCurrent.ML_N,...
            metadataCurrent.z_resolution_nm,metadataCurrent.x_resolution_nm, metadataCurrent.FourierOrders,metadataCurrent.plSelect);
        outputdata=grazing_angle_deg;
        outputdata=ef.inc_top_reflected.efficiency;
        outputdata=90-ef.inc_top_reflected.theta;
                
        idx_DesignOrder=find(metadataCurrent.FourierOrders:-1:0==metadataCurrent.GR_order);
        eff=ef.inc_top_reflected.efficiency(idx_DesignOrder);% diffraction order wants to display
        
        theta=grazing_angle_deg;
                elapsedTime = toc;
        disp(['Simulation time at step ' num2str(i) ': ' datestr(seconds(elapsedTime), 'HH:MM:SS')]);
        
    
else
while  k>0 && loopIndices{1} <= loopBounds(1)
    % Your code inside the innermost loop
    disp(['Loop indices: ' num2str(cell2mat(loopIndices))]);
    
    metadataCurrent=metadata;
    
    % transfer current value to metadataCurrent.
    for idx_i=1:loopNumber
        loopIndices_i=loopIndices{idx_i};
        metadataCurrent.(fields{loop.Idx{idx_i}})=metadata.(fields{loop.Idx{idx_i}})(loopIndices_i);
        outputdata1(idx_outputdata0,idx_i)=metadataCurrent.(fields{loop.Idx{idx_i}});
    end
%     if metadataCurrent.ML_d_nm==5 && metadataCurrent.photonEnergy_eV>=7500
    metadataCurrent.grBA_deg=0.15042*metadataCurrent.ML_d_nm + 0.00475;
    metadataCurrent.grAntiBA_deg=metadataCurrent.grBA_deg*3;
    
    thetaEst=estimateTheta(metadataCurrent.material_HZ,metadataCurrent.material_LZ, ...
        metadataCurrent.ML_order, metadataCurrent.ML_d_nm,metadataCurrent.ML_d_HZtod, ...
        metadataCurrent.GR_order,metadataCurrent.GR_groove,metadataCurrent.grPeriod_lpermm,metadataCurrent.photonEnergy_eV);

    i=1;
    for grazing_angle_deg=thetaEst %(thetaEst-0.06*2:0.005:thetaEst+0.06*2)
        ef = efficiency_bgrML(metadataCurrent.grPeriod_lpermm, ...
            metadataCurrent.grBA_deg, metadataCurrent.grAntiBA_deg, metadataCurrent.photonEnergy_eV, grazing_angle_deg, ...
            metadataCurrent.material_sub,metadataCurrent.material_HZ,metadataCurrent.material_LZ, ...
            metadataCurrent.ML_d_nm,metadataCurrent.ML_d_HZtod,metadataCurrent.ML_N,...
            metadataCurrent.z_resolution_nm,metadataCurrent.x_resolution_nm, metadataCurrent.FourierOrders,metadataCurrent.plSelect);
        outputdata(i,1)=grazing_angle_deg;
        outputdata(i,((1:length(ef.inc_top_reflected.efficiency)))*2)=ef.inc_top_reflected.efficiency;
        outputdata(i,((1:length(ef.inc_top_reflected.efficiency)))*2+1)=90-ef.inc_top_reflected.theta;
                
        idx_DesignOrder=find(metadataCurrent.FourierOrders:-1:0==metadataCurrent.GR_order);
        eff(i)=ef.inc_top_reflected.efficiency(idx_DesignOrder);% diffraction order wants to display
        
        theta(i)=grazing_angle_deg;
                elapsedTime = toc;
        disp(['Simulation time at step ' num2str(i) ': ' datestr(seconds(elapsedTime), 'HH:MM:SS')]);
        i=i+1;        
    end
    plot(theta,eff);xlabel('theta');title('Diffraction efficiency');legend('TE');ylabel('-1th diffraction efficiency');pause(eps);            
    hold on    
    metadataCurrent.outputRaw=outputdata;
    
    idx_MaxEff=find(outputdata(:,idx_DesignOrder*2)==max(outputdata(:,idx_DesignOrder*2)));
    metadataCurrent.Ginc_deg=outputdata(idx_MaxEff,1);
    metadataCurrent.Diff_deg=outputdata(idx_MaxEff,idx_DesignOrder*2+1);
    metadataCurrent.MaxEff_GrOrder=outputdata(idx_MaxEff,idx_DesignOrder*2);
    outputdata0{idx_outputdata0}=metadataCurrent;
    outputdata1(idx_outputdata0,idx_i+1:idx_i+3)=[metadataCurrent.Ginc_deg,metadataCurrent.Diff_deg,metadataCurrent.MaxEff_GrOrder]; 
    % Stop the timer at the end of the simulation
    elapsedTime = toc;
    disp(['Total simulation time: ' datestr(seconds(elapsedTime), 'HH:MM:SS')]);
    
    idx_outputdata0=idx_outputdata0+1;
%     end
    % Update loop indices
    k = numLoops;
    while k > 0
        loopIndices{k} = loopIndices{k} + 1;
        if loopIndices{k} <= loopBounds(k)
            break;
        else
            loopIndices{k} = 1;
            k = k - 1;
        end
    end
end
end
%     save(['bgrCmp_effAll_v6_04',num2str(metadata.grPeriod_lpermm),'LinePerMm.mat'],'outputdata0');
%     save(['bgrCmp_effAll_v6_04_2',num2str(metadata.grPeriod_lpermm),'LinePerMm.mat'],'outputdata1');



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