addpath(genpath('../'))
load('Images/Masks','Masks');
load('Images/Regions','Regions');
load('Images/Backgrounds','Backgrounds');

%% Importation des gold standards + var bool�ennes
load('Gold_Standards.mat');
Import=true;
Log_scale=true;

Mask2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3;
Im2Train=[1]; 

%% Initialisation des param�tres pour le test de S�rensen-Dice
pow_min=-5; pow_max=5; pow_step=1;
Nbins=[2,5,10,15,20];

pow_length=length(pow_min:pow_step:pow_max);
bins_length=length(Nbins);

Lambda=nan(1,pow_length);
Dice_Histo=nan(3,bins_length,pow_length);
Times=nan(3,bins_length,pow_length);

%% Param�tres num�riques
itermax=200; 
mu=0.1; beta=0.5; theta=0;
stop_1=1.e-6; stop_2=1.e-6;

% Arguments utilis�s pour voir les histogrammes 
fig=figure('Name','Histogramme des r�gions � segmenter','NumberTitle','off');
hist_visibility='off';

%% Lecture des images
Image1=zeros(248);
Image1(62:186,62:186)=255;
Image2=double(imread('eight.tif'));
Image3=double(imread('TumeurCerveaubis.png'));

NbImages=length(Im2Train);  %Nombre d'images test�es
TumeurBis=3; %Num�ro de l'image contenant le d�tail de la tumeure
Images={Image1,Image2,Image3};
%% Normalisation des images
for im=Im2Train
    Images{im}=Image_Normalisation(Images{im},"2D");
end

% Importation des masques ou cr�ation de ces derniers
if ~Import
    for m=Mask2Change
        fprintf("Choisissez la r�gion � segmenter sur l'image %d\n",m);
        Masks{m}=roipoly(Images{m});
    end
    close;
    
    for b=Background2Change
        fprintf("Choisissez un �l�ment de fond de l'image %d\n",b);
        Backgrounds{b}=roipoly(Images{b});
    end
    close;
    
    for r=Region2Change
        fprintf("Choisissez la r�gion de la radio � comperer avec le gold standard\n");
        Regions{r}=roipoly(Images{r});
    end
    close;
end

fprintf('%d it�rations pr�vues\n',pow_length);
iter=0;

%% Boucle des lambdas

figure(fig);
ax1=subplot(2,1,1); ax2=subplot(2,1,2);
for ind=pow_min:pow_step:pow_max
    if Log_scale
        lambda=10^ind;
    else
        lambda=ind;
    end
    iter=iter+1;
    fprintf("D�but de l'it�ration %d\n",iter);
    
    Lambda(iter)=lambda;
    
    %% Boucles des images + nombres de bins dans l'histogramme
    for im=Im2Train
        for nb=1:bins_length
            [g0,g1,T,sigma_1,sigma_2,sigma_3,b]=create_histo(Images{im},Masks{im},Backgrounds{im},Nbins(nb),ax1,ax2,hist_visibility);
            
            tic
            [Ub_Histo,~, ~,~, ~]=histo_loco(double(Masks{im}),g0,g1,b,T,sigma_1,sigma_2,sigma_3,lambda,mu,beta,theta,stop_1,stop_2,itermax);
            Times(im,nb,iter)=toc;
            
            if im==TumeurBis
                Ub_Histo=Ub_Histo & Regions{TumeurBis};
            end
            Dice_Histo(im,nb,iter)=dice(Ub_Histo,Gold_Standards{im});
        end
        fprintf("fin de l'it�ration %d.%d\n",iter,im);
    end
end
close(fig);

%% Affichage des r�sultats des tests de S�rensen-Dice
f1=figure('Name','Resultats des tests de S�rensen-Dice','NumberTitle','off');
f2=figure('Name','Etude du temps de calcul','NumberTitle','off');
f3=figure('Name','Images','NumberTitle','off');

%% Cr�ation de la l�gende 
leg={};
for nb=1:bins_length
    leg{nb}=['n',num2str(nb),' = ',num2str(Nbins(nb)),' bins'];
end
leg=reshape(leg,bins_length,1);

iter=0;
for im=Im2Train
    iter=iter+1;
    figure(f1);
    
    subplot(1,NbImages,iter);
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
        title({"Evolution de l'indice de Dice fonction du param�tre \lambda";['Image ',num2str(im)]});
    else
        title(["Image ",num2str(im)]);
    end
    
    figure(f2)
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
    if im==ceil(NbImages/2)
        title({'Evolution du temps de calcul en fonction du param�tre \lambda';['Image ',num2str(im)]});
    else
        title(["Image ",num2str(im)]);
    end
    
    figure(f3)
    subplot(1,NbImages,iter)
    imshow(Images{im});
    hold on
    contour(Masks{im},'r','Linewidth',3);
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