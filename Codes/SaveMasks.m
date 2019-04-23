function SaveMasks(filename, Foregrounds,Foreground2Change,Backgrounds, Background2Change,Regions,Region2Change)
%% Fonction utilisee pour sauvegarder les masques modifies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Sauvegarde des masques modifies')
for f=Foreground2Change
    FG=Foregrounds{f};
    save(['Images/Foregrounds/',filename{f},'_FG.mat'],'FG');
end

for b=Background2Change
    BG=Backgrounds{b};
    save(['Images/Backgrounds/',filename{b},'_BG.mat'],'BG');
end

for r=Region2Change
    R=Regions{r};
    save(['Images/Regions/',filename{r},'_Region.mat'],'R');
end

end