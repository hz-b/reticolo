%   conique    exemple_V9_conique_anisotrope


clear;
LD=.8;% longueur d'onde
D=1.5;% pas du reseau
teta0=30;
nh=1.9;nb=1.433+2i;ro=nh*sin(teta0*pi/180);
delta0=40;
parm=res0;parm.not_io=1; 
nn=20;

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures=cell(1,4);
textures{1}= nh;
% anisotropie generale non amplificatrice:  il faut que -1i*(epsilon-epsilon') soit definie >0
epsilon=randn(3)+1i*randn(3);epsilon=epsilon+epsilon';H=randn(3,1)+1i*randn(3,1);epsilon=1i*H*H'+epsilon';
nprov=rand;
parm.res1.change_index={ [nprov,  .1+5i,1.7+2i,1.3+1i] };% anisotropie diagonale


textures{2}=nb;
textures{3}={ epsilon} ;
textures{4}={ [-D/4,D/4],[nprov,1]} ;

% initialisation
[aa,neff]=res1(LD,D,textures,nn,ro,delta0,parm);
% visualisation du  mode 14 de la texture 4
res1(aa,neff, 4, 14);


% definition du profil et calcul de la diffraction
profil={[.2,.5,.7,.2], [1,3,4,2]};

ef=res2(aa,profil);

% reflexion
R_TE=sum(ef.TEinc_top_reflected.efficiency)
R_TM=sum(ef.TMinc_top_reflected.efficiency)


x=linspace(-D(1)/2,D(1)/2,51);
einc= ef.TEinc_top.PlaneWave_TE_Eu;
parm.res3.trace=1; ;%trace automatique


[e,z,o]=res3(x,aa,profil,einc,parm);
% Pertes
    
parm.res3.pertes_poynting=1;


x=[-D(1)/2,D(1)/2];parm.res3.trace=0;
[e,z,o,w,PP,P,p,XX,wx,Flux_Poynting]=res3(x,aa,profil,einc,parm);
bilan_energie= 1-(sum(ef.TEinc_top_reflected.efficiency)+sum(ef.TEinc_top_transmitted.efficiency))-(sum(PP)-Flux_Poynting(end))/(.5*prod(D))
% en tenant compte des pertes dans le substrat -Flux_Poynting(end)/(.5*prod(D))




