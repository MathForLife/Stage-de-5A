function Image=Image_Normalisation(Image,Dimension)
[~,~,p]=size(Image);

if strcmp(Dimension,"2D") && p>1
    Image=Image(:,:,1);
end

Min=min(min(Image));
Max=max(max(Image));
    
Image=Image-Min;
Image=Image/(Max-Min);
end