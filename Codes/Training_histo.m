clear all; close all;

image_names={'BrainHole','BrainMeta_A','BrainTumor_A','BrainTumorDetail','Coins',...
    'Flag','GBM_A','GeometricShape','Gliome003_S','Lung','Parenchyme_C','pneumopath_6_A','Square'};extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Train=[12]; NbImages=length(Im2Train);
ImWithRegion=7:8; 

import_masks=true; change_masks=false;  Bruitage=false; CompareTexture=true; import_textures=true; change_textures=false; choose_texture=true;
Foreground2Change=[6]; Background2Change=[6]; Region2Change=[]; Texture2Change=[];

[Images, Textures, Foregrounds, Backgrounds,Regions,Gold_Standards]=ImportImageMasks(image_names,extension,Im2Train,import_masks,ImWithRegion,change_masks,Foreground2Change,Background2Change,Region2Change,import_textures,choose_texture,change_textures,Texture2Change);
%% Initialisation des parametres du test de Sorensen-Dice
fprintf('Initialisation des parametres \n');
pow_min=-5; pow_max=5; pow_step=1;
Log_scale=true; 

NbinsF=[2,5,10]; NbinsB=[2,5,10];

pow_length=length(pow_min:pow_step:pow_max);
bins_length=length(NbinsF); 

Lambda=nan(1,pow_length);
Dice_Histo=nan(NbImages,bins_length,pow_length);

%% Cas ou on compare l'image BrainTumorDetail avec ses 2 Gold Standards
if ismember('BrainTumorDetail',image_names(Im2Train))
    TestGSBrainTumor=true;
    load(['Gold_Standards/',image_names{6},'_GS2.mat'],'GS');
    GS_FullTumor=GS;
    Dice_FullTumor=nan(bins_length,pow_length);
else
    TestGSBrainTumor=false;
end
Times=nan(NbImages,bins_length,pow_length);

%% Parametres numeriques (se referer au programme HrstogrammeSegmentation pour le detail de ces parametres)
itermax=500;  
mu=0.5; beta=0.5; theta=1; epsilon=1;
stop_u=-1.e-8; stop_J=-1.e-8;

% Arguments utilises pour voir les histogrammes 
fig=figure('Name','Histogramme des regions a segmenter','NumberTitle','off');
hist_visibility='off'; cumulative.value=true;
if cumulative.value
    cumulative.normalisation='cdf';
else
    cumulative.normalisation='probability';
end

fprintf('%d iterations prevues\n',pow_length);
iter=0;


%% Creation des figures 
figure(fig);
ax1=subplot(2,1,1); ax2=subplot(2,1,2);

PlotOptions={ax1,ax2,hist_visibility};
StopConditions=[itermax,stop_u,stop_J];
Parameters=[mu, beta, theta, epsilon];

%% Boucle des lambdas
for ind=pow_min:pow_step:pow_max
    if Log_scale
        lambda=10^ind;
    else
        lambda=ind;
    end
    iter=iter+1;
    fprintf("Debut de l'iteration %d\n",iter);
    
    Lambda(iter)=lambda;
    
    %% Boucles des images + nombres de bins dans l'histogramme
    for im=Im2Train
        for bins=1:bins_length     
            Nbins=[NbinsF(bins),NbinsB(bins)];
            tic
            [Ub_Histo,~, ~,~,~,~,~,~]=HistogrammeSegmentation(Images{im},lambda,Foregrounds{im},Backgrounds{im},Nbins,cumulative,hist_visibility,Parameters,StopConditions,Textures{im});
            Times(im,bins,iter)=toc;
            
            if ismember(im,ImWithRegion)
                Ub_Histo=Ub_Histo & Regions{im};
            end
%             figure(2);
%             subplot(141); imshow(Gold_Standards{im});
%             subplot(142); imshow(Ub_Histo_temp);
%             subplot(143); imshow(Ub_Histo);
%             subplot(144); imshow(Regions{TumeurBis});
            
            Dice_Histo(im,bins,iter)=SorensenDice(Ub_Histo,Gold_Standards{im});
            if strcmp(image_names{im},'BrainTumorDetail')
                Dice_FullTumor(bins,iter)=SorensenDice(Ub_Histo,GS_FullTumor);
            end
        end
        fprintf("fin de l'iteration %d.%d\n",iter,im);
    end
end
close(fig);

%% Affichage des resultats des tests de Sorensen-Dice
f1=figure('Name','Resultats des tests de Sorensen-Dice','NumberTitle','off');
f2=figure('Name','Evolution du temps de calcul','NumberTitle','off');
f3=figure('Name','Images','NumberTitle','off');

%% Creation de la legende 
leg={}; 
for bins=1:bins_length
    leg{bins}=['Nf = ',num2str(NbinsF(bins)),' bins, Nb = ',num2str(NbinsB(bins)),' bins'];
end
leg=leg'; 

if TestGSBrainTumor
    NbDiceplot=NbImages+1;
else
    NbDiceplot=NbImages;
end

iter=0; iterDiceplot=0;
for im=Im2Train
    iter=iter+1; iterDiceplot=iterDiceplot+1;
    
    figure(f1);
    subplot(1,NbDiceplot,iterDiceplot);
    Dice=reshape(Dice_Histo(im,:,:),bins_length,pow_length);
    
    if Log_scale
        semilogx(Lambda,Dice);
    else
        plot(Lambda,Dice);
    end
    
    legend(leg,'Location','best');
    xlabel('\lambda');
    ylabel('Dice index');
    title(image_names{im});
    
    if strcmp(image_names{im},'BrainTumorDetail') 
        iterDiceplot=iterDiceplot+1;
        subplot(1,NbDiceplot,iterDiceplot);
        
        if Log_scale
            semilogx(Lambda,Dice_FullTumor);
        else
            plot(Lambda,Dice_FullTumor);
        end
        
        legend(leg,'Location','best');
        xlabel('\lambda');
        ylabel('Dice index');
        title([image_names{im},'\_GS2']);
    end
    figure(f2);
    subplot(1,NbImages,iter);
    t=reshape(Times(im,:,:),bins_length,pow_length);
    if Log_scale
        semilogx(Lambda,t);
    else
        plot(Lambda,t);
    end
    legend(leg,'Location','best');
    
    ylabel('t');
    xlabel('\lambda');
    title(image_names{im});
    
    figure(f3);
    subplot(1,NbImages,iter)
    imshow(Images{im});
    hold on
    contour(Foregrounds{im},'r','Linewidth',3);
    contour(Backgrounds{im},'g','Linewidth',3);
    hold off
     if im==ceil(NbImages/2)
        title({'Image + Initial mask + Background mask';image_names{im}});
        legend('Initial segmentation mask','Background mask');
        legend('Initial \partial \Omega_1','Initial \partial \Omega_0');
    else
        title(image_names{im});
     end
end