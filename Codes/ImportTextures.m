function Textures=ImportTextures(image_names,Im2Select,texture_names,Text2Select)
fprintf('Importation des textures');
Textures=struct([]);
for im=Im2Select
    for text=Text2Select
        name=texture_names{text};
        
        load(['../Images/Textures/',name,'/',image_names{im},'_TX.mat'],name);
        switch name
            case 'Energy'
                Textures(im).Energy=Energy;
            case 'Entropy'
                Textures(im).Entropy=Entropy;
            case 'Correlation'
                Textures(im).Correlation=Correlation;
            case 'IDM'
                Textures(im).IDM=IDM;
            case 'Inertia'
                Textures(im).Inertia=Inertia;
            case 'Cluster_Shade'
                Textures(im).Cluster_Shade=Cluster_Shade;
            case 'Cluster_Prominence'
                Textures(im).Cluster_Prominence=Cluster_Prominence;
        end
    end
end