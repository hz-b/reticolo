%   2  D    exemple_V9_0D_anisotrope


clear;
LD=.5;% longueur d'onde
D=1.5;% pas du reseau n'importe quoi
teta0=30;
nh=1.9;ro=nh*sind(teta0);
delta0=0;
parm=res0;parm.not_io=1; 
nn=0; % si pas de struture 

% description des textures y compris le substrat et le superstrat et les milieux homogenes
textures=cell(1,3);
textures{1}= nh;  
epsilon=[[2.1160         0    0.7165];[0    1.3995         0]; [0.7165         0    2.1160]];
textures{2}={ epsilon} ;
textures{3}=1.433;
% initialisation
[aa,neff]=res1(LD,D,textures,nn,ro,delta0,parm);
% definition du profil et calcul de la diffraction
h=.5;
profil={[.2,h,.2], [1,2,3]};

ef=res2(aa,profil);

% reflexion
R_TE=ef.TEinc_top_reflected.efficiency
R_TM=ef.TMinc_top_reflected.efficiency


x=linspace(-D(1)/2,D(1)/2,51);y=0;% x,y coordonnes du maillage en x et y (les cordonnes en z sont determinees par le calcul)
einc=[1,0];   %  composantes du champ e incident dans le repere u v  par defaut le champ incident vient du haut
parm=res0;parm.not_io=1;  %parametres par defaut
parm.res3.trace=1 ;%trace automatique


[e,z,o]=res3(x,y,aa,profil,einc,parm);
