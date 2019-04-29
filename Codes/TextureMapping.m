function Texture=TextureMapping(Image,Image_name,d_patch,d_glcm,ChooseTexture)
fprintf('Creation des cartes de texture\n');
sz=size(Image);

Energy=zeros(sz); Entropy=zeros(sz); Correlation=zeros(sz); IDM=zeros(sz);
Inertia=zeros(sz); Cluster_Shade=zeros(sz); Cluster_Prominence=zeros(sz);
Texture_names={'Color','Energy','Entropy','Correlation','IDM','Inertia','Cluster_Shade','Cluster_Prominence'};

[mat_C,mat_L]=meshgrid(1:8);
offsets=[0 , d_glcm ; d_glcm , d_glcm ; d_glcm , 0 ; d_glcm , -d_glcm];
%% Boucles sur les pixels de l'image
for i=1:sz(1)
    % 
    i_min=max(i-d_patch,1); i_max=min(i+d_patch,sz(1));
    for j=1:sz(2)        
        j_min=max(j-d_patch,1); j_max=min(j+d_patch,sz(2));
        Patch=Image(i_min:i_max,j_min:j_max);
        
        GLCM=graycomatrix(Patch,'Offset',offsets,'Symmetric',true);

        for k=1:4
            glcm_temp=GLCM(:,:,k)/((i_max-i_min)*(j_max-j_min));
            mask_temp=glcm_temp~=0;
            
            mu=sum(sum(mat_L.*glcm_temp));
            var=sum(sum(glcm_temp.*(mat_L-mu).^2));
            
            Energy(i,j)=Energy(i,j)+sum(sum(glcm_temp.^2));
            
            Entropy(i,j)=Entropy(i,j)+sum(sum(glcm_temp(mask_temp).*log2(glcm_temp(mask_temp))));
            
            Correlation(i,j)=Correlation(i,j)+sum(sum((mat_L-mu).*glcm_temp.*(mat_C-mu)))/var;
            
            IDM(i,j)=IDM(i,j)+sum(sum(glcm_temp./(1+(mat_L-mat_C).^2)));
            
            Inertia(i,j)=Inertia(i,j)+sum(sum((mat_L-mat_C).^2.*glcm_temp));
            
            Cluster_Shade(i,j)=Cluster_Shade(i,j)+sum(sum(((mat_L-mu)+(mat_C-mu)).^3.*glcm_temp));
            
            Cluster_Prominence(i,j)=Cluster_Prominence(i,j)+sum(sum(((mat_L-mu)+(mat_C-mu)).^4.*glcm_temp));
        end
%         Energy(i,j)=Energy(i,j)/4;
%         Entropy(i,j)=Entropy(i,j)/4;        
%         Correlation(i,j)=Correlation(i,j)/4;
%         IDM(i,j)=IDM(i,j)/4;
%         Inertia(i,j)=Inertia(i,j)/4;        
%         Cluster_Shade(i,j)=Cluster_Shade(i,j)/4;
%         Cluster_Prominence(i,j)=Cluster_Prominence(i,j)/4;
        
    end
end
%% Renormalisation des indices de texture moyennes
Energy=Image_Normalisation(Energy/4,"2D");
Entropy=Image_Normalisation(Entropy/4,"2D");
Correlation=Image_Normalisation(Correlation/4,"2D");
IDM=Image_Normalisation(IDM/4,"2D");
Inertia=Image_Normalisation(Inertia/4,"2D");
Cluster_Shade=Image_Normalisation(Cluster_Shade/4,"2D");
Cluster_Prominence=Image_Normalisation(Cluster_Prominence/4,"2D");
        
%Cluster_Prominence=adapthisteq(Cluster_Prominence);
%% Affichage des differents indices de texture 
figure(20)
ax1=subplot(331); imshow(Image); title('Image');
ax2=subplot(332); imagesc(Energy); axis off; title('Energy'); colorbar();
ax3=subplot(333); imagesc(Entropy); axis off; title('Entropy'); colorbar();
ax4=subplot(334); imagesc(Correlation); axis off; title('Correlation'); colorbar();
ax5=subplot(335); imagesc(IDM); axis off; title('IDM'); colorbar();
ax6=subplot(336); imagesc(Inertia); axis off; title('Inertia'); colorbar();
ax7=subplot(337); imagesc(Cluster_Shade); axis off; title('Cluster Shade'); colorbar();
ax8=subplot(338); imagesc(Cluster_Prominence); axis off; title('Cluster Prominence'); colorbar();
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'xy');

%% Choix des textures a prendre en compte 
if ChooseTexture
    Text2Consider=input(['Rentrer le nom des indicateurs de texture a prendre en compte sous forme d un vecteur \n',...
        '1 : Color, 2 : Energy, 3 : Entropy, 4 : Correlation\n',...
        '5 : IDM, 6 : Inertia, 7 : Cluster_Shade, 8 : Cluster_Prominence\n']);
else
    Text2Consider=1:8;
end

%% Creation de la structure Texture comportant les differents indicateurs de texture
for text=Text2Consider
    name=Texture_names{text};
    if ~strcmp(name,'Color')
        save(['../Images/Textures/',name,'/',Image_name,'_TX.mat'],name);
    end
    
    switch name
        case 'Color'
            Texture.Color=Image;
        case 'Energy'
            Texture.Energy=Energy;
        case 'Entropy'
            Texture.Entropy=Entropy;
        case 'Correlation'
            Texture.Correlation=Correlation;
        case 'IDM'
            Texture.IDM=IDM;
        case 'Inertia'
            Texture.Inertia=Inertia;
        case 'Cluster_Shade'
            Texture.Cluster_Shade=Cluster_Shade;
        case 'Cluster_Prominence'
            Texture.Cluster_Prominence=Cluster_Prominence;
    end
end
