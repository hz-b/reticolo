function y=retinv_cell(x)
%  function y=retinv_cell(x);
% inversion 'analytique' de cell_array carre (jusque à n=3)
n=size(x,1);y=x;
switch(n);
	
case 1;y{1}=1./x{1};
	
case 2;prv=x{1,1}.*x{2,2}-x{1,2}.*x{2,1};
y{1,1}=x{2,2}./prv;y{2,2}=x{1,1}./prv;y{1,2}=-x{1,2}./prv;y{2,1}=-x{2,1}./prv;	

case 3;prv=x{1,1}.*x{2,2}.*x{3,3}+x{2,1}.*x{3,2}.*x{1,3}+x{1,2}.*x{2,3}.*x{3,1}...
-x{3,1}.*x{2,2}.*x{1,3}-x{2,1}.*x{1,2}.*x{3,3}-x{3,2}.*x{2,3}.*x{1,1};		

y{1,1}=(x{2,2}.*x{3,3}-x{3,2}.*x{2,3})./prv;	
y{1,2}=(x{1,3}.*x{3,2}-x{1,2}.*x{3,3})./prv;	
y{1,3}=(x{1,2}.*x{2,3}-x{2,2}.*x{1,3})./prv;

y{2,1}=(x{2,3}.*x{3,1}-x{2,1}.*x{3,3})./prv;	
y{2,2}=(x{1,1}.*x{3,3}-x{3,1}.*x{1,3})./prv;	
y{2,3}=(x{1,3}.*x{2,1}-x{1,1}.*x{2,3})./prv;

y{3,1}=(x{2,1}.*x{3,2}-x{2,2}.*x{3,1})./prv;	
y{3,2}=(x{1,2}.*x{3,1}-x{1,1}.*x{3,2})./prv;	
y{3,3}=(x{1,1}.*x{2,2}-x{2,1}.*x{1,2})./prv;	

otherwise;
prv=zeros(n);
for k=1:length(x{1,1}(:));
for ii=1:n;for jj=1:n;prv(ii,jj)=x{ii,jj}(k);end;end;
prv=inv(prv);
for ii=1:n;for jj=1:n;y{ii,jj}(k)=prv(ii,jj);end;end;
end;

end