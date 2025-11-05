%BMG_RCWA
clear;retio

%以nm为单位
%膜层厚度读取，计算平均膜厚
E_keV=2.5;
wavelength=1.2398./E_keV;
N= 50 ;  %50 达 1mm
NN=N*2;
d_spacing=4;
gamma=0.4;

% 
%profile from AFM
Prf=importdata('1723_PEAXIS_2400_G7_ghosts_pos9_00054-2D.txt');
Prf=Prf*1E9; % m to nm
% for i=1:1
%     A=[];
% A(:,1)=[Prf(:,1);Prf(:,1)+max(Prf(:,1))+0.001];
% A(:,2)=[Prf(:,2);Prf(:,2)];
% Prf=A;
% end
% %%% or generate ideal one
% line=2400;%l/mm
% LD=1000./line.*1000;
% 
% BA=1*pi./180;%弧度
% antiBA=3*BA;%弧度

% GD=LD.*sin(BA).*sin(antiBA)./sin(BA+antiBA);
% D_BA=GD./tan(BA);
% Prf=[0,0;D_BA,GD;LD,0];

% generate grid %unit nm
x_L=Prf(end,1); x_step=10; x=linspace(0,x_L,round(x_L/x_step)+1);%
z_L=max(Prf(:,2))+d_spacing*N+0.1; z_step=0.05; z=linspace(z_L,0,round(z_L/z_step)+1);%;
LD=Prf(end,1);
[X,Z]=meshgrid(x,z);
% interplot the z of profilr based on x
Prf_z=interp1(Prf(:,1)',Prf(:,2)',x);
%%%折射率轮廓%%%%%
%材料折射率参数读取
n_inc=1;
n_real_A=1-0.00021810481;n_imag_A=1.9496438E-05;
n_real_S=1-7.4526305E-05;n_imag_S=1.3561066E-06;
n_real_sub=1-7.5815149E-05;n_imag_sub=1.4569904E-05;
n_A=n_real_A+n_imag_A*1i;
n_S=n_real_S+n_imag_S*1i;
n_sub=n_real_sub+n_imag_sub*1i;
n=Z.*0;
Zz=n; % to store the layer profile
% generate model, add Optical constant
for i=1:N
    Prf0=Prf_z+d_spacing*(i-1);
    Prf1=Prf_z+d_spacing*(i-1)+d_spacing*gamma;
    Prf2=Prf_z+d_spacing*i;
    if i==1
        P=find(Z<Prf0);         Zz(P)=1; n(P)=n_sub;%substrate
        P=find(Z>=Prf0&Z<Prf1);  Zz(P)=2; n(P)=n_A; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  Zz(P)=3; n(P)=n_S; % spacing layer
    elseif i==N %top layer
        P=find(Z>=Prf2);         Zz(P)=4; n(P)=n_inc;% background
        P=find(Z>=Prf0&Z<Prf1);  Zz(P)=2; n(P)=n_A; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  Zz(P)=3; n(P)=n_S; % spacing layer     
    else
        P=find(Z>=Prf0&Z<Prf1);  Zz(P)=2; n(P)=n_A; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  Zz(P)=3; n(P)=n_S; % spacing layer     
    end
end

imagesc(Zz)
% colormap(gca, 'jet')
colorbar;
% rowdist=ones(1,round(z_L/z_step)+1);
rowdist=ones(1,round(z_L/z_step)+1);
textures_tmp= [mat2cell(X,rowdist) mat2cell(n,rowdist)];
for i = 1:round(z_L/z_step)+1
    textures{i} = {textures_tmp{i,1}, textures_tmp{i,2}};
end
% texture_list=1:round(z_L/z_step)+1;
% thicknessstep_list=ones(1,round(z_L/z_step)+1).*z_step;
texture_list=1:round(z_L/z_step)+1;
thicknessstep_list=ones(1,round(z_L/z_step)+1).*z_step;
profile={thicknessstep_list,texture_list};%每层0.1nm，层序号top开始为1


%RCWA
te=[];tm=[];theta=[];
plotprofile=1;
for theta_0=linspace(87.9,88.3,4)
    theta_0
    nn=5;% ordres de fourier
    k_parallel=n_inc*sin(theta_0*pi/180);
    for te_tm=[1,-1]
        parm=res0(te_tm);  %res0(1):TE;res0(-1):TM;% initialisation des parametres par defaut
        aa=res1(wavelength,LD,textures,nn,k_parallel,parm);
        % aa = res1(wavelength,period,textures,nn,k_parallel,parm)
        
        %plotprofile
%         if plotprofile==1
%             figure;
%             x=linspace(-LD,LD,2*num_X+1);
%             [e,z,index]=res3(x,aa,profile,1,parm);
%             %[e,z,index] = res3(x,aa,profile,inc,parm)% Computation of the electromagnetic fields%profile1
%             retcolor(x,z,real(index));xlabel('X');ylabel('Z');title('profile');axis equal;pause(eps)
%             plotprofile=0;
%         end
        
        % -1级次衍射效率
        result=res2(aa,profile);%result = res2(aa, profile)
        if te_tm==1
            te=[te,result.inc_top_reflected.efficiency{-1}];
            theta=[theta,theta_0];
        else
            tm=[tm,result.inc_top_reflected.efficiency{-1}];
            avg_tetm=(te+tm)./2;
            plot(theta,te,theta,tm,'--',theta,avg_tetm,'*');xlabel('theta');title('Diffraction efficiency');legend('TE','TM','AVG');ylabel('-1th diffraction efficiency');pause(eps);
        end
    end
end

retio

%%%得到结果并读取参数及结构%%%%
%%%优化算法程序%%%%
%%%将新参数进行处理回到第一步%%%%

