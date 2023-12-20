%    1 D    exemple_V9_1D_anisotrope


clear;
LD=.8;% longueur d'onde
D=1.5;% pas du reseau
teta0=30;
nh=1.9;nb=1.433;ro=nh*sin(teta0*pi/180);
for pol=[1,-1];
parm=res0(pol);parm.not_io=1; 
nn=100;

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures=cell(1,3);
textures{1}= nh;
nprov=rand;
parm.res1.change_index={ [nprov,  .1+5i,1.7+2i,1.3+1i] };

textures{2}=nb;
textures{3}={ [-D/4,D/4],[nprov,1.5]} ;

% initialisation
[aa,neff]=res1(LD,D,textures,nn,ro,parm);
% visualisation des 2 premiers modes de la texture 4
res1(aa,neff, 3, 1);res1(aa,neff, 3, 2);


% definition du profil et calcul de la diffraction
profil={[.2,.5,.2], [1,3,2]};

ef=res2(aa,profil);
% reflexion
R=sum(ef.inc_top_reflected.efficiency)


x=linspace(-D(1)/2,D(1)/2,51);
if pol==1;einc=ef.inc_top.PlaneWave_E(2);else;einc=ef.inc_top.PlaneWave_H(2);end;
parm.res3.trace=1; %trace automatique
parm.res3.npts=[[10,10,10];[1,4,1]];  % points en z
parm.res3.gauss_x=30;% points en x
[e,z,o]=res3(x,aa,profil,einc,parm);
% Pertes en calculant l'integrale
    


x=[-D(1)/2,D(1)/2];parm.res3.trace=0;
[e,z,o,w,PP,P,p,XX,wx]=res3(x,aa,profil,einc,parm);

bilan_energie= 1-(sum(ef.inc_top_reflected.efficiency)+sum(ef.inc_top_transmitted.efficiency))-sum(PP)/(.5*D)


end;% pol

