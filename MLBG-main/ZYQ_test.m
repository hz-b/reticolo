clear;
%BMG_RCWA
clear;retio

%以nm为单位
%膜层厚度读取，计算平均膜厚
E_keV=2.5;
wavelength=1.2398./E_keV;
N=100;
NN=N*2;
d_spacing=5.8;
gamma=0.4;

% 
% %profile from AFM
% Prf=importdata('20230912_Hefei_1DprofileST1.txt');
% Prf=Prf*1000; % um to nm
%%% or generate ideal one
line=2400;%l/mm
LD=1000./line.*1000;

BA=1*pi./180;%弧度
antiBA=3*BA;%弧度

GD=LD.*sin(BA).*sin(antiBA)./sin(BA+antiBA);
D_BA=GD./tan(BA);
Prf=[0,0;D_BA,GD;LD,0];

% generate grid %unit nm
x_L=Prf(end,1); x_step=5; x=linspace(0,x_L,round(x_L/x_step)+1);%
z_L=max(Prf(:,2))+d_spacing*N+1; z_step=0.05; z=linspace(z_L,0,round(z_L/z_step)+1);%;
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
    Prf0=Prf_z+1+d_spacing*(i-1);
    Prf1=Prf_z+1+d_spacing*(i-1)+d_spacing*gamma;
    Prf2=Prf_z+1+d_spacing*i;
    if i==1
        P=find(Z<Prf0);         Zz(P)=4; n(P)=n_sub;%substrate
        P=find(Z>=Prf0&Z<Prf1);  Zz(P)=2; n(P)=n_A; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  Zz(P)=3; n(P)=n_S; % spacing layer
    elseif i==N %top layer
        P=find(Z>=Prf2);         Zz(P)=1; n(P)=n_inc;% background
        P=find(Z>=Prf0&Z<Prf1);  Zz(P)=2; n(P)=n_A; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  Zz(P)=3; n(P)=n_S; % spacing layer     
    else
        P=find(Z>=Prf0&Z<Prf1);  Zz(P)=2; n(P)=n_A; % absorption layer
        P=find(Z>=Prf1&Z<Prf2);  Zz(P)=3; n(P)=n_S; % spacing layer     
    end
end

imagesc(Zz)


%%if mod(XX,2)==1 %奇数处理 else %偶数处理
%第1层为最顶层，第NN层为最底层
for j=1:NN
    if mod(j,2)==1
        d(j)=d_spacing*(1-gamma);%共NN层，奇数层S厚度
    else
        d(j)=d_spacing*gamma;%偶数层A厚度
    end
end

%匹配闪耀角（P2、P3max），单独开一个函数m文件
line=2400;%l/mm
LD=1000./line.*1000;

L_total=sum(d)+GD;

%%%折射率轮廓%%%%%
%材料折射率参数读取

%生成层数据
thicknessstep=0.05;%每层厚度0.1nm
num_Z=round(z_L/z_step);%分层数num_Z,
L_sum=num_Z.*thicknessstep;%整个结构高度
num_X=round(x_L/x_step);
X_step=LD./num_X;
X=linspace(0,LD,num_X+1);
Z=linspace(0,L_sum,num_Z+1);

%各点材料判定公式
for xj=1:num_X+1
    if (xj*X_step)<=D_BA
        deltaz(xj)=(D_BA-X(xj)).*tan(BA)+thicknessstep;
    else
        deltaz(xj)=(X(xj)-D_BA).*tan(antiBA)+thicknessstep;
    end
    for zj=1:num_Z+1
        if Z(zj)<deltaz(xj)
            n(zj,xj)=n_inc; n_map(zj,xj)=1;
        else
            n(zj,xj)=n_S;n_map(zj,xj)=3;
        end
        for Nj=2:NN
            if Z(zj)>=(sum(d(1:Nj-1))+deltaz(xj))&&Z(zj)<(sum(d(1:Nj))+deltaz(xj))
                if mod(Nj,2)==1%共NN层，奇数层S
                    n(zj,xj)=n_S;n_map(zj,xj)=3;
                else
                    n(zj,xj)=n_A;n_map(zj,xj)=2;%偶数层A
                end
            end
        end
        if Z(zj)>=sum(d)+deltaz(xj)
            n(zj,xj)=n_sub;n_map(zj,xj)=4;
        end
    end
end

imagesc(n_map);
colormap jet;
colorbar;

Z=1:100;
Y=1:30;
dimt=ones(30,1)';
L=meshgrid(Z,Y);
LL=L.*3;
T=mat2cell(L,dimt);
TT=mat2cell(LL,dimt);
x=[1,2];
y=1:30;
pl=meshgrid(x,y);

cell21= @(n) {L(n,:),LL(n,:)};
x = cell21(1:5)
