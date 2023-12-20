function z=retprod_cell(varargin)
% function y=retprod_cell(x1,x2,x3,..);
% produit de cell_array 'matrices'

if nargin==1;z=varargin{1};return;end;
if nargin>2;z=retprod_cell(retprod_cell(varargin{1:2}),varargin{3:end});return;end;
n=size(varargin{1},1);[m,p]=size(varargin{2});z=cell(n,p);[z{:}]=deal(0);[z{:}]=deal(zeros(size(varargin{1}{1})));
for in=1:n;for ip=1:p;for im=1:m;z{in,ip}=z{in,ip}+varargin{1}{in,im}.*varargin{2}{im,ip};end;end;end;
