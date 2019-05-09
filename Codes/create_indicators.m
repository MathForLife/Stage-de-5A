function [I,c1,c0]=create_indicators(Texture, Foreground, Background, vect_length)


Texture_name=fieldnames(Texture); NbTexture=length(Texture_name);
I=zeros(NbTexture,vect_length); 
c1=zeros(NbTexture,1); c0=zeros(NbTexture,1); 

for text=1:NbTexture
    switch Texture_name{text}
        case 'Color'
            IndicateurTexture=Texture.Color;
        case 'Energy'
            IndicateurTexture=Texture.Energy;
        case 'Entropy'
            IndicateurTexture=Texture.Entropy;
        case 'Correlation'
            IndicateurTexture=Texture.Correlation;
        case 'IDM'
            IndicateurTexture=Texture.IDM;
        case 'Inertia'
            IndicateurTexture=Texture.Inertia;
        case 'Cluster_Shade'
            IndicateurTexture=Texture.Cluster_Shade;
        case 'Cluster_Prominence'
            IndicateurTexture=Texture.Cluster_Prominence;
    end
    I(text,:)=reshape(IndicateurTexture,1,vect_length);
    c1(text)=mean(IndicateurTexture(Foreground));
    c0(text)=mean(IndicateurTexture(Background));
end