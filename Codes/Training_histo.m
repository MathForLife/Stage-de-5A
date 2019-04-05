clear all; close all;

filename={'Square','GeometricShape','Coins','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Train=[1,3,5]; 
ImWithRegion=4:7; 

ChangeMasks=false;
Foreground2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3;

[Images, Foregrounds, Backgrounds,Regions,Gold_Standards]=ImportImageMasks(filename,extension,Im2Train,ImWithRegion,ChangeMasks,Foreground2Change,Background2Change,Region2Change);
% Cas où on compare l'image BrainTumorDetail avec ses 2 GS
if ismember(5,Im2Train)
    load(['Gold_Standards/',filename{5},'_GS2.mat'],'GS');
    GS_FullTumor=GS;
end
%% Initialisation des parametres du test de Sorensen-Dice
fprintf('Initialisation des parametres \n');
pow_min=-6; pow_max=6; pow_step=1;
Log_scale=true;

NbinsF=[20]; NbinsB=[10];

pow_length=length(pow_min:pow_step:pow_max);
bins_length=length(NbinsF); NbImages=length(Im2Train);

Lambda=nan(1,pow_length);
Dice_Histo=nan(NbImages,bins_length,pow_length);
Dice_FullTumor=nan(bins_length,pow_length);
Times=nan(NbImages,bins_length,pow_length);
%% Parametres numeriques (se referer au programme HrstogrammeSegmentation pour le detail de ces parametres)
itermax=200;  
mu=0.5; beta=0.5; theta=1; epsilon=10;
stop_u=-1.e-8; stop_J=-1.e-8;

% Arguments utilises pour voir les histogrammes 
fig=figure('Name','Histogramme des regions a segmenter','NumberTitle','off');
hist_visibility='off'; cumulative=false;

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
            [Ub_Histo,~,~,~,~]=HistogrammeSegmentation(Images{im},Foregrounds{im},Backgrounds{im},Nbins,cumulative,lambda,Parameters,PlotOptions,StopConditions);
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
            if strcmp(filename{im},'BrainTumorDetail')
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
leg={}; leg_tumor={};
for bins=1:bins_length
    leg{bins}=['Nf = ',num2str(NbinsF(bins)),' bins, Nb = ',num2str(NbinsB(bins)),' bins'];
    leg_tumor{bins}=['Pf = ',num2str(NbinsF(bins)),' bins, Pb = ',num2str(NbinsB(bins)),' bins'];
end
leg=leg'; leg_tumor=leg_tumor';
iter=0;
for im=Im2Train
    iter=iter+1;
    
    figure(f1);
    subplot(1,NbImages,iter);
    Dice=reshape(Dice_Histo(im,:,:),bins_length,pow_length);
    if strcmp(filename{im},'BrainTumorDetail') 
        Dice=[Dice;Dice_FullTumor];
    end
    
    if Log_scale
        semilogx(Lambda,Dice);
    else
        plot(Lambda,Dice);
    end
    
    if strcmp(filename{im},'BrainTumorDetail') 
        legend([leg;leg_tumor],'Location','best');
    else
        legend(leg,'Location','best');
    end
    xlabel('\lambda');
    ylabel('Dice index');
    title(filename{im});
    
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
    title(filename{im});
    
    figure(f3);
    subplot(1,NbImages,iter)
    imshow(Images{im});
    hold on
    contour(Foregrounds{im},'r','Linewidth',3);
    contour(Backgrounds{im},'g','Linewidth',3);
    hold off
     if im==ceil(NbImages/2)
        title({'Image + Initial mask + Background mask';['Image ',num2str(im)]});
        legend('Initial segmentation mask','Background mask');
        legend('Initial \partial \Omega_1','Initial \partial \Omega_0');
    else
        title(["Image ",num2str(im)]);
    end
    title(["Image et masques initiaux ",num2str(im)]);
end

%% Sauvegarde des masques modifies
if ChangeMasks
    SaveMasks(filename, Foregrounds,Foreground2Change,Backgrounds, Background2Change,Regions,Region2Change)
end