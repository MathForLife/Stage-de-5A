function M=div(varargin)
%calcul de la divergence d'un vecteur de dimension 2 ou 3
%px et py de meme taille
%(de sorte que div=-(grad)^*)
%Syntaxe: div(px,py)
px=varargin{1};
py=varargin{2};
if nargin==3
    pz=varargin{3};
end

[m,n,p]=size(px);

M=zeros(m,n,p);
Mx=M; My=M; Mz=M;

Mx(2:m-1,:,:)=px(2:m-1,:,:)-px(1:m-2,:,:);
Mx(1,:,:)=px(1,:,:);
Mx(m,:,:)=-px(m-1,:,:);

My(:,2:n-1,:)=py(:,2:n-1,:)-py(:,1:n-2,:);
My(:,1,:)=py(:,1,:);
My(:,n,:)=-py(:,n-1,:);

if nargin==3
    Mz(:,:,2:n-1)=pz(:,:,2:p-1)-pz(:,:,1:p-2);
    Mz(:,:,1)=pz(:,:,1);
    Mz(:,:,n)=-pz(:,:,p-1);
end
M=Mx+My+Mz;



