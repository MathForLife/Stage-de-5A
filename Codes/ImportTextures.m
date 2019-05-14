function Textures=ImportTextures(Images,image_names,Im2Consider,ChooseTexture)
fprintf('Importation des textures\n');
Texture_names={'Color','Energy','Entropy','Correlation','IDM','Inertia','Cluster_Shade','Cluster_Prominence'};
Textures={};

for im=Im2Consider
    Image=Images{im};
    
    for text=2:8
        name=Texture_names{text};
        load(['../Images/Textures/',name,'/',image_names{im},'_TX.mat'],name);
    end
    
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
        Text2Select=input(['Rentrer l indice des indicateurs de texture a prendre en compte sous forme d un vecteur \n',...
            '1 : Color, 2 : Energy, 3 : Entropy, 4 : Correlation\n',...
            '5 : IDM, 6 : Inertia, 7 : Cluster_Shade, 8 : Cluster_Prominence\n']);
    else
        Text2Select=1;
    end
    
    for text=Text2Select
        name=Texture_names{text};
        
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
    Textures{im}=Texture;    
end

end