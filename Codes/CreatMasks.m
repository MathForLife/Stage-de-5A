function [Foregrounds, Backgrounds,Regions]=CreatMasks(Images,image_names,extension,Foregrounds,Foreground2Change,Backgrounds,Background2Change,Regions,Region2Change)
%% Creation de nouveaux masques
fprintf('Creation de nouveaux masques\n');

for f=Foreground2Change
    fprintf("Choisissez la region a segmenter sur l'image %s\n",image_names{f});
    Foregrounds{f}=roipoly(Images{f});
    % Sauvegarde du nouveau masque
    imwrite(Foregrounds{f},['Images/Foregrounds/',image_names{f},'_FG',extension]);
end
close;

for b=Background2Change
    fprintf("Choisissez un element de fond de l'image %d\n",image_names{b});
    Backgrounds{b}=roipoly(Images{b});
    % Sauvegarde du nouveau masque
    imwrite(Backgrounds{b},['Images/Backgrounds/',image_names{b},'_BG',extension]);
end
close;

for r=Region2Change
    fprintf("Choisissez la region de la radio a comparer avec le gold standard\n");
    Regions{r}=roipoly(Images{r});
    % Sauvegarde du nouveau masque
    imwrite(Regions{r},['Images/Regions/',image_names{r},'_Region',extension]);
end
close;

end
