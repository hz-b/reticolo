%BMG_RCWA
% clear;
retio

%以nm为单位
%膜层厚度读取，计算平均膜厚
E_keV=3;
wavelength=1.2398./E_keV;
N= 60 ;  %50 达 1mm
NN=N*2;
d_spacing=6.5;
gamma=0.45;


%profile lamellar grating
%%% or generate ideal one
line=150;%l/mm
p_nm=1000./line.*1000;
grDepth=15;%nm
grRatioWidSp=0.65;
sub=10;%nm

%estimate Inc
ML_order=1;
GR_groove=1;
GR_order=1;
%RCWA
te=[];tm=[];theta=[];
plotprofile=1;
% for theta_0=min(thetaEst-0.12,89.35):0.005:min(thetaEst+0.1,89.95)
energy=380;
n_inc=1;

for i=1:length(energy)   
    energy(i)
    wavelength=1239.8/energy(i);
    [profile,textures,Prf]=buildModule(p_nm,grDepth,grRatioWidSp,sub,energy(i));
    theta_0=90-2;  
    k_parallel=n_inc*sind(theta_0);
    nn=5;% ordres de fourier
    for te_tm=[1,-1]
        parm=res0(1);  %res0(1):TE;res0(-1):TM;% initialisation des parametres par defaut
        aa=res1(wavelength,p_nm,textures,nn,k_parallel,parm);
        result=res2(aa,profile);%result = res2(aa, profile)
%                     show if needed
            x=linspace(-p_nm,p_nm,501);% on trace 2 periodes
            parm.res3.trace=1 ; % trace automatique
            parm.res3.cale=[];
            parm.res3.npts=[10,80,10];
            [e,z,o]=res3(x,aa,profile,1,parm);
            axis square
            set(gcf,'WindowStyle','docked')
        if te_tm==1
            te=[te,result.inc_top_reflected.efficiency{-GR_order*GR_groove}];
            theta=[theta,theta_0];
        else
            tm=[tm,result.inc_top_reflected.efficiency{-GR_order*GR_groove}];
            avg_tetm=(te+tm)./2;
        end
%         efficiencyCurve(i,1)=energy(i);
%         efficiencyCurve(i,2)=te;
%         efficiencyCurve(i,3)=tm;
%         efficiencyCurve(i,4)=avg_tetm;
    end

end
figure
plot(energy,te)   
hold on
plot(energy,tm) 
plot(energy,avg_tetm) 
retio

    figure
    plot(Prf(:,1),Prf(:,2))
%%%得到结果并读取参数及结构%%%%
%%%优化算法程序%%%%
%%%将新参数进行处理回到第一步%%%%

function [profile,textures,Prf]=buildModule(LD,grDepth,grRatioWidSp,sub,energy)
    Prf(:,1)=[0,0, LD*grRatioWidSp,LD*grRatioWidSp,LD,LD];
    Prf(:,2)=[-sub,grDepth, grDepth,0,0,-10];


    % generate grid %unit nm
    AA_speed=2;
    x_L=Prf(end,1); x_step=5*AA_speed^0.5; x=linspace(0,x_L,round(x_L/x_step)+1);%
    z_L=max(Prf(:,2))+5; z_step=0.1*AA_speed; z=linspace(z_L,-sub,round(z_L/z_step)+1);%;
    LD=Prf(end,1);
    [X,Z]=meshgrid(x,z);
    % interplot the z of profilr based on x
    idx=find(x<=LD*grRatioWidSp);
    Prf_z(idx)=grDepth;
    idx=find(x>LD*grRatioWidSp);
    Prf_z(idx)=0;
    %%%折射率轮廓%%%%%
    %材料折射率参数读取
    n_subFile=importdata('20231219_cxro_Au_Density19.32.txt');n_subFile=n_subFile.data;
    n_inc=1;
    n_real_sub=1-interp1(n_subFile(:,1),n_subFile(:,2),energy);n_imag_sub=interp1(n_subFile(:,1),n_subFile(:,3),energy);
    n_sub=n_real_sub+n_imag_sub*1i;
    n=Z.*0;

    Zz=n; % to store the layer profile
    % generate model, add Optical constant
    P=find(Z<Prf_z);         Zz(P)=1; n(P)=n_sub;%substrate
    %top layer
    P=find(Z>=Prf_z);         Zz(P)=4; n(P)=n_inc;% background

    % imagesc(Zz)
    % colormap(gca, 'jet')
    % colorbar;
    rowdist=ones(1,round(z_L/z_step)+1);
    textures_tmp= [mat2cell(X,rowdist) mat2cell(n,rowdist)];
    for i = 1:round(z_L/z_step)+1
        textures{i} = {textures_tmp{i,1}, textures_tmp{i,2}};
    end
    texture_list=1:round(z_L/z_step)+1;
    thicknessstep_list=ones(1,round(z_L/z_step)+1).*z_step;
    profile={thicknessstep_list,texture_list};%每层0.1nm，层序号top开始为1
end

