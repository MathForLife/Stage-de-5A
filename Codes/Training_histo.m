clear all; close all;

filename={'Square','GeometricShape','Coins','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
extension='.png'; addpath('../Images/');

%% Choix des images sur lesquelles entrainer les algos + importation des masques
Im2Train=[1,2,3]; Im2Reg=4:7; 
Images={}; Foregrounds={}; Backgrounds={}; Regions={}; Gold_Standards={};
disp('Importation des images et des masques')
for im=Im2Train
    Images{im}=imread([filename{im},extension]);
    load(['Foregrounds/',filename{im},'_FG.mat'],'FG');
    Foregrounds{im}=FG;
    load(['Backgrounds/',filename{im},'_BG.mat'],'BG');
    Backgrounds{im}=BG;
    if ismember(im,Im2Reg)
        load(['Regions/',filename{im},'_Region.mat'],'R');
        Regions{im}=R;
    end    
    load(['Gold_Standards/',filename{im},'_GS.mat'],'GS');
    Gold_Standards{im}=GS;
end

%% Possibilite de modifier les masque d'initialisation en indiquant les numero des images correspondant au masque
ChangeMasks=false;
Mask2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3;

%% Initialisation des parametres du test de Sorensen-Dice
disp('Initialisation des parametres ')
pow_min=1; pow_max=1; pow_step=1;
Log_scale=true;

NbinsF=[2,4]; NbinsB=[2,4];

pow_length=length(pow_min:pow_step:pow_max);
bins_length=length(NbinsF); nbIm=length(Im2Train);

Lambda=nan(1,pow_length);
Dice_Histo=nan(nbIm,bins_length,pow_length);
Times=nan(nbIm,bins_length,pow_length);
%% Parametres numeriques (se referer au programme HostogrammeSegmentation pourle detail de ces parametres)
itermax=200; 
mu=0.5; beta=0.5; theta=1; epsilon=1;
stop_u=-1.e-6; stop_J=-1.e-6;

% Arguments utilises pour voir les histogrammes 
fig=figure('Name','Histogramme des regions a segmenter','NumberTitle','off');
hist_visibility='on';

%% Normalisation des images
disp('Normalisation des images')
for im=Im2Train
    Images{im}=Image_Normalisation(Images{im},"2D");
end

%% Creation de nouveaux masques 
if ChangeMasks
    disp('Modification des masques initiaux')
    for f=Foreground2Change
        fprintf("Choisissez la region a segmenter sur l'image %d\n",f);
        Foregrounds{f}=roipoly(Images{f});
    end
    close;
    
    for b=Background2Change
        fprintf("Choisissez un element de fond de l'image %d\n",b);
        Backgrounds{b}=roipoly(Images{b});
    end
    close;
    
    for r=Region2Change
        fprintf("Choisissez la region de la radio a comparer avec le gold standard\n");
        Regions{r}=roipoly(Images{r});
    end
    close;
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
        for nb=1:bins_length     
            Nbins=[NbinsF(nb),NbinsB(nb)];
            tic
            [Ub_Histo,~,~,~,~]=HistogrammeSegmentation(Images{im},Foregrounds{im},Backgrounds{im},Nbins,lambda,Parameters,PlotOptions,StopConditions);
            Times(im,nb,iter)=toc;
            
            if ismember(im,Im2Reg)
                Ub_Histo=Ub_Histo & Regions{im};
            end
%             figure(2);
%             subplot(141); imshow(Gold_Standards{im});
%             subplot(142); imshow(Ub_Histo_temp);
%             subplot(143); imshow(Ub_Histo);
%             subplot(144); imshow(Regions{TumeurBis});
            
            Dice_Histo(im,nb,iter)=SorensenDice(Ub_Histo,Gold_Standards{im});
        end
        fprintf("fin de l'iteration %d.%d\n",iter,im);
    end
end
%close(fig);

%% Affichage des resultats des tests de Sorensen-Dice
f1=figure('Name','Resultats des tests de Sorensen-Dice','NumberTitle','off');
%f2=figure('Name','Etude du temps de calcul','NumberTitle','off');
f2=figure('Name','Images','NumberTitle','off');

%% Creation de la legende 
leg={};
for nb=1:bins_length
    leg{nb}=['n',num2str(nb),' = ',num2str(NbinsF(nb)),' bins'];
end
leg=reshape(leg,bins_length,1);

iter=0;
for im=Im2Train
    iter=iter+1;
    figure(f1);
    
    subplot(2,NbImages,iter);
    Dice=reshape(Dice_Histo(im,:,:),bins_length,pow_length);
    if Log_scale
        semilogx(Lambda,Dice);
    else
        plot(Lambda,Dice);
    end
    legend(leg,'Location','best');

    xlabel('\lambda');
    ylabel('Dice index');
    if im==ceil(NbImages/2)
        title({"Evolution de l'indice de Dice en fonction du parametre \lambda";['Image ',num2str(im)]});
    else
        title(["Image ",num2str(im)]);
    end
    
    subplot(2,NbImages,iter+NbImages);
    t=reshape(Times(im,:,:),bins_length,pow_length);
    if Log_scale
        semilogx(Lambda,t);
    else
        plot(Lambda,t);
    end
    legend(leg,'Location','best');
    
    ylabel('t');
    xlabel('\lambda');
    if im==ceil(NbImages/2)
        title({'Evolution du temps de calcul en fonction du parametre \lambda';['Image ',num2str(im)]});
    else
        title(["Image ",num2str(im)]);
    end
    
    figure(f2)
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
    title(["Image et masque initial ",num2str(im)]);
end


if ~Import
    save('Images/Masks','Masks');
    save('Images/Backgrounds','Backgrounds');
    save('Images/Regions','Regions');
end