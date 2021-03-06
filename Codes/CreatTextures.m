function Textures=CreatTextures(Images,Textures,image_names,Texture2Change,choose_texture)
fprintf('Creation de nouvelles textures\n');

for t=Texture2Change
    if choose_texture
        d_patch=ceil(input('Rentrer la valeur de la distance entre le centre du patch et une de ses extremite\n'));
        d_glcm=ceil(input('Rentrer la valeur de la norme des vecteurs utilises pour la calcul de la matrice de coocurence\n'));
    else
        d_patch=5; d_glcm=1;
    end
    Textures{t}=TextureMapping(Images{t},image_names{t},choose_texture,d_patch,d_glcm);
end