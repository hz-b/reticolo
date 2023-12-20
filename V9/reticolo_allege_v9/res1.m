function varargout=res1(varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               %%%%%%%%%%%%%%%%%
%               %     2 D       %
%               %%%%%%%%%%%%%%%%%
%    function [a,nef]=res1(LD,D,TEXTURES,nn,ro,delta0,parm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% construction de a  cell array contenant des donnees sur les textures
% et nef cell array des 'indices efficaces' des textures (pour les metaux infiniment conducteurs nef{ii}=[])
%  ces 'indices efficaces' sont classes par attenuation decroissante
%
%  LD,D  longueur d'onde ,pas en unites metriques
%  TEXTURES ={ { n1,[cx1,cy1,dx1,dy1,ni1,k1],[cx2,cy2,dx2,dy2,ni2 k2],...[cxn,cyn,dxn,dyn,nin kn]}  % premiere texture
%      ,   {n2,...   }  % deuxieme texture
%      ,   {nm,...    }  } % derniere texture
%
%  TEXTURES{2} ={ n1, [cx1,cy1,dx1,dy1,ni1,k1],[cx2,cy2,dx2,dy2,ni2,k2],...[cxn,cyn,dxn,dyn,nin,kn]}
%    n1:indice de la base
%   [cx1,cy1,          dx1,dy1,             ni1     k1   ]: premiere inclusion
%    centre    largeurs en x et y      indice
%  k1=1  l'inclusion est un rectangle  de  cotes  dx1,dy1
%  k1 >1  l'inclusion est une ellipse de grands axes  dx1,dy1,  approchee par k1 rectangles 
%                   la surface totale de l'inclusion etant celle de l'ellipse
%     si le rectangle ou l'ellipse a une dimension plus grande que le pas il y a chevauchement (indice ni1 dans la partie commune)
%
%
%  on peut aussi au lieu d'une base homogene avoir un reseau  1 D
%  decrit par un tableau de points de discontinuites suivi du tableau des indices a gauche des points
%  les points de discontinuitee doivent etre en ordre croissant sur un intervalle strictement inferieur au pas 
%   et etre au moins 2  
%      si ces tableaux sont des vecteurs ligne le reseau est invariant en y (cas conique)
%      si ces tableaux sont des vecteurs colonnes le reseau est invariant en x 
%
%   ATTENTION: l'ordre des inclusions est important car elles s'ecrasent l'une l'autre ..
% 
%   on peut aussi definir des plaques de metal infiniment conducteur percees de trous rectangulaires NE SE CHEVAUCHANT PAS
%  TEXTURES{ }= { inf, [cx1,cy1,dx1,dy1,ni1],[cx2,cy2,dx2,dy2,ni2],..} pour le metal electrique
%  TEXTURES{ }= {-inf, [cx1,cy1,dx1,dy1,ni1],[cx2,cy2,dx2,dy2,ni2],..} pour le metal magnetique
%   [cx1,cy1,          dx1,dy1,         ni1       ]: premier trou
%    centre    largeurs en x et y      indice
%  par exemple  TEXTURES{ }= { inf} est le metal massif en haut ou en bas 
% 
% 
%  Possibilite d'introduire des couches homogenes (cad non structurées) avev un epsilon 3 X 3 ou meme un mu
%  TEXTURE{4}={epsilon}  ou TEXTURE{4}={epsilon,mu} epsilon mu matrices 3X3
%  
%
% SIMPLIFICATIONS D'ECRITURE
%  si on a une seule texture il n'est pas necessaire de la mettre en tableau de cell array
%   TEXTURES ={ n1, [cx1,cy1,dx1,dy1,ni1,k1],[cx2,cy2,dx2,dy2,ni2 k2],...[cxn,cyn,dxn,dyn,nin kn]}
%  pour les milieux homogenes on peut entrer  TEXTURES{k} =n1  
%
%
%
%  nn nombre de termes de fourier en x et y
%           si size(nn)=[1,2]    -nn(1) a nn(1) en x -nn(2) a nn(2) en y
%           si size(nn)=[2,2]     nn(1,1) a nn(2,1) en x  nn(1,2) a nn(2,2) en y
%           si size(nn)=[1,1]    -nn(1) a nn(1) en x  0  en y (cas conique)
%           si size(nn)=[2,1]    nn(1,1) a nn(2,1) en x  0  en y (cas conique) 
%
%  ro,deltao: incidence beta0=beta0=ro*[cos(delta0*pi/180),sin(delta0*pi/180)]   par defaut ro=1, delta0=0 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  symetries
%  parm.sym.x=x0   x=x0 :plan de symetrie seulement si beta0(1)=0     par defaut parm.sym.x=[]  pas de symetrie en x 
%  parm.sym.y=y0   y=y0 :plan de symetrie seulement si beta0(2)=0     par defaut parm.sym.y=[]  pas de symetrie en y
%  parm.sym.pol=1;   1 TE  -1:TM  par defaut 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  parm.res1.trace:  1  trace des textures     (defaut 0)
%
%  parm.res1.nx    parm.res1.xlimite (valeurs de x et nombre de points pour le trace)
%  parm.res1.ny    parm.res1.ylimite (valeurs de y et nombre de points pour le trace)
%     (par defaut la maille du reseau centree avec 100*100 points)
%  exemple:  parm.res1.xlimite=[-D(1),D(1)]; parm.res1.nx=500; % parametres pour le trace des textures en x
%
%  parm.res1.angles: pour les ellipses   1 :regulier en angles  0 meilleure repartition reguliers en surface (defaut 1)
%  parm.res1.calcul: 1   calcul  (defaut 1)  (si on ne veut que visualiser les textures prendre parm.res1.calcul=0 )
%  parm.res1.champ:  options pour le trace futur des champs  1: les champs sont bien calcules( gourmand en temps et memoire)
%        0 :les champs non continus sont calcules de facon approchee  et Ez et Hz ne sont pas calcules et remplaces par 0
%                              (defaut 0)
%  parm.res1.ftemp: 1   fichiers temporaires  (defaut 1) 
%  parm.res1.fperm: 'aa'   le resultat est mis sur un fichier permanent  de nom 'aa123..'  (defaut [] donc pas ecriture0) 
%
% si neff est en sortie neff est un cell array de la dimension de TEXTURES qui contient les indices effectifs de chague texture

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               %%%%%%%%%%%%%%%%%
%               %     1 D       %
%               %%%%%%%%%%%%%%%%%
%    function [a,nef]=res1D(LD,D,TEXTURES,nn,pol,beta0,parm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% construction de a  cell array contenant des donnees sur les textures  cas 1 D
% et nef cell array des 'indices efficaces' des textures (pour les metaux infiniment conducteurs nef{ii}=[])
%  ces 'indices efficaces' sont classes par attenuation decroissante
%
%  LD,D  longueur d'onde ,pas en unites metriques
%  TEXTURES ={ {x1,n1},{x2,n2},...{xn,nn}} 
%
%
%  profil des tranches du reseau
%  decrites par un tableau de points de discontinuites suivi du tableau des indices a gauche des points
%  les points de discontinuitee doivent etre en ordre croissant sur un intervalle strictement inferieur au pas 
%   et etre au moins 2  
%
% 
%   on peut aussi definir des tranches de metal infiniment conducteur percees de trous rectangulaires NE SE CHEVAUCHANT PAS
%  TEXTURES{ }= { inf, [cx1,dx1,ni1],[cx2,dx2,ni2],..} pour le metal electrique
%  TEXTURES{ }= {-inf, [cx1,dx1,n1],[cx2,dx2,ni2],..} pour le metal magnetique
%   [cx1,          dx1,         ni1       ]: premier trou
%    centre    largeurs      indice
%  par exemple  TEXTURES{ }= { inf} est le metal massif en haut ou en bas  
%
% SIMPLIFICATIONS D'ECRITURE
%  si on a une seule texture il n'est pas necessaire de la mettre en tableau de cell array

%  pour les milieux homogenes on peut entrer  TEXTURES{k} =n1  
%
%
%
%  nn nombre de termes de fourier
%
%           si length(nn)=1    -nn a nn 
%           si length(nn)=2    nn(1) a nn(2)
%
%   beta0=n*sin(teta*pi/180)  par defaut beta0=0 
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  symetries
%  parm.sym.x=x0   x=x0 :plan de symetrie    par defaut parm.sym.x=[]  pas de symetrie 
%     possible uniquement si beta0=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polarisation  pol 1:TE  -1:TM  par defaut TE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%   AIGUILLAGE   %
%%%%%%%%%%%%%%%%%%

if iscell(varargin{1});[varargout{1:nargout}]=trace_mode(varargin{:});return;end;

if nargin>=6 & isstruct(varargin{6});[varargout{1:nargout}]=res1_1D(varargin{:});% cas 1 D
else;
[varargout{1:nargout}]=res1_2D(varargin{:});% cas 2 D
end;    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,nef]=res1_2D(LD,D,TEXTURES,nn,ro,delta0,parm);
if nargin<7;parm=[];end;if isempty(parm);parm=res0;end;
if parm.not_io==1;[prv,vmax]=retio([],inf*i);parm.res1.ftemp=0;end;
    
if nargin<5;ro=0;end;
if nargin<6;delta0=0;end;beta0=ro*[cos(delta0*pi/180),sin(delta0*pi/180)];
UN=LD/(2*pi);
if (length(D)<2);D=[D,D];parm.sym.y=0;end;if size(nn,2)==1;nn=[nn,zeros(size(nn))];end; % cas conique( symetrie en y)
if size(nn,1)==1;nn=[-nn;nn];end; % en general les termes de fourier vont de -nn a nn mais on se reserve la possibilitee generale

sym=[0,0,0,0]; % symetries
if abs(sin(delta0*pi/180))<10*eps;sym(1:2)=-parm.sym.pol;end;
if abs(cos(delta0*pi/180))<10*eps;sym(1:2)=parm.sym.pol;end;
if ~isempty(parm.sym.x)&(abs(beta0(1))<10*eps)&(nn(2,1)==-nn(1,1));sym(3)=parm.sym.x/UN;else;sym(1)=0;end;   % symetrie par rapport au plan x=x0 
if ~isempty(parm.sym.y)&(abs(beta0(2))<10*eps)&(nn(2,2)==-nn(1,2));sym(4)=parm.sym.y/UN;else;sym(2)=0;end;   % symetrie par rapport au plan y=y0 

if parm.res1.trace==1;    % trace 
if isempty(parm.res1.xlimite);x=linspace(-D(1)/2,D(1)/2,parm.res1.nx)/UN;else;x=linspace(parm.res1.xlimite(1),parm.res1.xlimite(2),parm.res1.nx)/UN;end;    
if isempty(parm.res1.ylimite);y=linspace(-D(2)/2,D(2)/2,parm.res1.ny)/UN;else;y=linspace(parm.res1.ylimite(1),parm.res1.ylimite(2),parm.res1.ny)/UN;end;    
figure;
end;

% mise en forme de TEXTURE
if ~iscell(TEXTURES);TEXTURES={TEXTURES};end;%                              TEXTURES = 1.5

if ~any(cellfun('isclass',TEXTURES,'cell'))&~all(cellfun('length',TEXTURES)==1);TEXTURES={TEXTURES};end;
%       TEXTURES={1.5,[..]}  ou    TEXTURES={[-.5,1,3],[2,1.3,1.5] , [0,0,6,6,  1,  5] ,  [0,0,2,2,  1.5,  1]    };
%   mais pas TEXTURES={1, 1.5 } qui signifie 2 milieux homogenes


NTEXTURES=size(TEXTURES,2);
nf=ceil(sqrt(NTEXTURES)); % pour le trace
a=cell(NTEXTURES,1);nef=cell(NTEXTURES,1);
sog=parm.res1.sog;
anisotrope=zeros(1,NTEXTURES);
for ii=1:NTEXTURES;if ~iscell(TEXTURES{ii});TEXTURES{ii}={TEXTURES{ii}};end;
if ~isfinite(TEXTURES{ii}{1}(1,1));sog=0;end;
end;

[init,n]=retinit(D/UN,[nn(1,1)+1i*(1-sog),nn(2,1),nn(1,2),nn(2,2)],[beta0,delta0*pi/180],sym);init{end}.granet=0;

for ii=1:NTEXTURES; %construction des a
minc=2;
if (length(TEXTURES{ii})>1)&(length(TEXTURES{ii}{1})==length(TEXTURES{ii}{2}));  % la base est une texture 1D 
minc=3; 
end;
for in=minc:size(TEXTURES{ii},2);
if parm.res1.angles==1&length(TEXTURES{ii}{in})>5;TEXTURES{ii}{in}(6)=-TEXTURES{ii}{in}(6);end;
end;

if all(size(TEXTURES{ii}{1})==[3,3])% <<<< couche anisotrope homogene
anisotrope(ii)=1;
if length(TEXTURES{ii})==1;% epsilon seul
u={1,1,reshape([retcolonne(eye(3));retcolonne(TEXTURES{ii}{1})],1,1,18)};% on multiplie pas par k0
else;% epsilon et mu
u={1,1,reshape([retcolonne(TEXTURES{ii}{2});retcolonne(TEXTURES{ii}{1})],1,1,18)};% on multiplie pas par k0
end
else;                               % <<<<<
u=retu(init,[TEXTURES{ii},{1,UN}]); % k0=1 UN mise a l'echelle spatial
% installation eventuelle de l'anisotropie diagonale
for indice=1:length(parm.res1.change_index);
ep=ret2ep(parm.res1.change_index{indice}(1));% indice bidon
if length(parm.res1.change_index{indice})==4;% epsilon seul
epnv=[1;1;1;...
parm.res1.change_index{indice}(2)^2;parm.res1.change_index{indice}(3)^2;parm.res1.change_index{indice}(4)^2];	
else% epsilon et mu
epnv=[parm.res1.change_index{indice}(5)^2;parm.res1.change_index{indice}(6)^2;parm.res1.change_index{indice}(7)^2; ...
parm.res1.change_index{indice}(2)^2;parm.res1.change_index{indice}(3)^2;parm.res1.change_index{indice}(4)^2];	
end;
u=installe_anisotropie(u,ep,epnv);
end;
end;                               % <<<<<

if length(u)==3 & length(u{1})==1 & length(u{2})==1;champ=1;else;champ=parm.res1.champ;end; % pour les milieux homogenes
if parm.res1.trace==1; % trace
[ee,xy]=rettestobjet(init,u,-1,[],{x,y});
subplot(nf,nf,ii);retcolor(xy{1}*UN,xy{2}*UN,real(ee.'));axis equal;
title(['texture  ',num2str(ii)],'fontsize',7);xlabel('X','fontsize',7);ylabel('Y','fontsize',7);pause(eps);
end;
if parm.res1.calcul==1;       %  diagonalisation
if anisotrope(ii);
[a{ii},neff]=retcouche_bis(init,{u,struct('teta',0,'phi',0)},1);
else;
[a{ii},neff]=retcouche(init,u,champ+1i*parm.res1.ftemp,parm.res1.li);
end;
if parm.res1.ftemp==1;a{ii}=retio(a{ii},1);end;
if nargout>1;nef{ii}=i*neff;  % indices effectifs
[prv,iii]=sort(imag(nef{ii}));nef{ii}=nef{ii}(iii);
f=find(abs(imag(nef{ii}))<100*eps);[prv,iii]=sort(real(nef{ii}(f)));nef{ii}(f)=nef{ii}(f(iii));
end;
end;
end;  % boucle sur ii
a={a,init,n,UN,D,beta0,sym};
if ~isempty(parm.res1.fperm);a=retio(a,parm.res1.fperm,0);end; % ecriture eventuelle sur fichier permanent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,nef]=res1_1D(LD,D,TEXTURES,nn,beta0,parm);

pol=1-parm.dim;
UN=LD/(2*pi);
if length(nn)==1;nn=[-nn,nn];end; 

if abs(beta0)<10*eps&~isempty(parm.sym.x);sym=[1,parm.sym.x];else sym=[];end;

if parm.res1.trace==1;    % trace 
if isempty(parm.res1.xlimite);x=linspace(-D/2,D/2,parm.res1.nx)/UN;else;x=linspace(parm.res1.xlimite(1),parm.res1.xlimite(2),parm.res1.nx)/UN;end;    
figure;
end;


% mise en forme de TEXTURE
if ~iscell(TEXTURES);TEXTURES={TEXTURES};end;%                              TEXTURES = 1.5
if ~any(cellfun('isclass',TEXTURES,'cell'))&~all(cellfun('length',TEXTURES)==1);TEXTURES={TEXTURES};end;
%       TEXTURES={1.5,[..]}  
%       TEXTURES={1, 1.5 }  signifie 2 milieux homogenes


NTEXTURES=size(TEXTURES,2);
nf=ceil(sqrt(NTEXTURES)); % pour le trace


a=cell(NTEXTURES,1);nef=cell(NTEXTURES,1);
sog=parm.res1.sog;
for ii=1:NTEXTURES;if ~iscell(TEXTURES{ii});TEXTURES{ii}={TEXTURES{ii}};end;
if ~isfinite(TEXTURES{ii}{1}(1));sog=0;end;
end;

[init,n]=retinit(D/UN,[nn(1)+i*(1-sog),nn(2)],beta0,sym);init{end}.granet=0;

for ii=1:NTEXTURES; % construction des a
u=retu(init,[TEXTURES{ii},{pol,1,UN}]);
% installation eventuelle de l'anisotropie diagonale
for indice=1:length(parm.res1.change_index);
% change_index=[nbidon,nx,ny,nz,mx,my,mz]
ep=retep(parm.res1.change_index{indice}(1),pol);% indice bidon
if pol==0;% TE
%	epnv = [eps_z; mu_z; 1/mu_x];
if length(parm.res1.change_index{indice})==4;% epsilon seul
epnv=[parm.res1.change_index{indice}(3)^2;1;1];
else% epsilon et mu
epnv=[parm.res1.change_index{indice}(3)^2;parm.res1.change_index{indice}(7)^2;1/parm.res1.change_index{indice}(5)^2];
end;
else;      % TM
if length(parm.res1.change_index{indice})==4;% epsilon seul
%	epnv = [mu_y;  eps_z;  1/eps_x];
epnv=[1;parm.res1.change_index{indice}(4)^2;1/parm.res1.change_index{indice}(2)^2];
else% epsilon et mu
epnv=[parm.res1.change_index{indice}(6)^2;parm.res1.change_index{indice}(4)^2;1/parm.res1.change_index{indice}(2)^2];
end;
end;  % TE TM

u=installe_anisotropie(u,ep,epnv);
end;
if parm.res1.trace==1; % trace
[ee,x]=rettestobjet(init,u,-1,[],x);
subplot(nf,nf,ii);plot(x*UN,real(ee));axis([min(x*UN),max(x*UN),0,max(abs(real(ee)))*1.1]);
title(['texture  ',num2str(ii)],'fontsize',7);xlabel('X','fontsize',7);ylabel('indice','fontsize',7);pause(eps);
end;

if length(u{1})==1;champ=1;else;champ=parm.res1.champ;end; % pour les milieux homogenes

[a{ii},neff]=retcouche(init,u,champ);
if nargout>1;nef{ii}=i*neff;  % indices effectifs
[prv,iii]=sort(imag(nef{ii}));nef{ii}=nef{ii}(iii);
f=find(abs(imag(nef{ii}))<100*eps);[prv,iii]=sort(real(nef{ii}(f)));nef{ii}(f)=nef{ii}(f(iii));
end;
end;  % boucle sur ii    
a={a,init,n,UN,D,beta0,sym,pol};
if ~isempty(parm.res1.fperm);a=retio(a,parm.res1.fperm,0);end; % ecriture eventuelle sur fichier permanent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u=installe_anisotropie(u,ep,epnv);
if length(u)==3;             % 2D 
prv=reshape(u{3},[],6);	
f=find(sum(abs(prv-repmat(ep.',size(prv,1),1)),2)<100*eps);
prv(f,:)=repmat(retcolonne(epnv,1),length(f),1);
u{3}=reshape(prv,size(u{3}));
else;                               % 1D
f=find(sum(abs(u{2}-repmat(ep,1,size(u{2},2))))<100*eps);
u{2}(:,f)=repmat(epnv,1,length(f));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e,o,x,y]=trace_mode(aa,neff,Num_texture,Num_mode,x,y);
% Num_texture=2;Num_mode=28;
init=aa{2};UN=aa{4};a=retio(aa{1}{Num_texture});
if init{end}.dim==2;% 2D ou conique
if nargin<6;y=0;end
if nargin<5;
period=aa{5};x=linspace(-period(1)/2,period(1)/2,150);y=linspace(-period(2)/2,period(2)/2,152);
end
[prv,num]=min(abs(a{5}*1i-neff{Num_texture}(Num_mode)));% find the mode
sh=retb(init,a,1.e-6,0,num,[]);
sb=retb(init,a,-1.e-6,0,[],[]);
[e,z,wz,o]=retchamp(init,{a},sh,sb,1,{x/UN,y/UN},[0,1,1]);% compute the field
if nargout==0;rettchamp(e,o,x,y,z,[1:6,12,1i]);end;

else;               % 1D
period=aa{5};
if nargin<5;
x=linspace(-period/2,period/2,150);
end
[prv,num]=min(abs(a{5}*1i-neff{Num_texture}(Num_mode)));% find the mode
sh=retb(init,a,1.e-6,0,num,[]);
sb=retb(init,a,-1.e-6,0,[],[]);
[e,z,wz,o]=retchamp(init,{a},sh,sb,1,x/UN,[period,1,1]);% compute the field
pol=aa{end};
if nargout==0;
ee=permute(e,[2,1,3]);ee=ee(:,:,[1,3,2]);oo=permute(o,[2,1]);% pour tracer avec rettchamp
rettchamp(ee,oo,z,x,pol,[1:3,12,1i]);
end;
end;	