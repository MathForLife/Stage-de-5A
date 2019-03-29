clear all; close all;

addpath(genpath('../'))
load('Images/Masks','Masks');
load('Images/Regions','Regions');

%% Importation des gold standards + var bool�ennes
load('Gold_Standards.mat');
Import=true;
Log_scale=false;

Mask2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3;
Im2Train=[2]; 
%% Initialisation des param�tres pour le test de S�rensen-Dice
pow_min=1; pow_max=10; pow_step=1;

pow_length=round(abs(pow_max-pow_min)/pow_step+1);
Lambda=nan(1,pow_length);

Dice_CV=nan(3,pow_length);
Dice_CP=nan(3,pow_length);
Dice_CEN=nan(3,pow_length);
Time=nan(3,3,pow_length);

%% Param�tres num�riques
itermax=200; cinconnu=false;
mu=0.5; eps=1;
theta=1; beta=0.5;
order=2; eta=1; n=10; % paramètres num. de CV
tho_u=0.5; tho_z=0.25;
stop_1=-1.e-6; stop_2=-1.e-6;
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
for i=Im2Train
    Images{i}=Image_Normalisation(Images{i},"2D");
end

% Importation des masques ou cr�ation de ces derniers
if ~Import
    for m=Mask2Change
        fprintf("Choisissez la r�gion � segmenter sur l'image %d\n",m);
        Masks{m}=roipoly(Images{m});
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
for ind=pow_min:pow_step:pow_max
    if Log_scale
        lambda=10^ind;
    else
        lambda=ind;
    end
    iter=iter+1;
    fprintf("D�but de l'it�ration %d\n",iter);
    
    Lambda(iter)=lambda;
    %% Boucle sur les images
    for i=Im2Train
        c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(~Masks{i}));
        
        tic
        [Ub_CV, J_CV, err_u,err_J, Niter_CV]=ChanVese(Images{i}, double(Masks{i}), lambda, eps, eta,stop_1,stop_2,c1,c2, cinconnu, n, order, itermax);
        t1=toc;
        
        tic
        [Ub_CP, ~, ~,~, ~]=DualFormulation(Images{i}, double(Masks{i}), lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
        t3=toc;
        
        tic
        [Ub_CEN, ~, ~,~, ~]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
        t2=toc;
        
        Time(i,:,iter)=[t1,t2,t3];
        
        if i==TumeurBis
            Ub_CV=Ub_CV & Regions{TumeurBis};
            Ub_CP=Ub_CP & Regions{TumeurBis};
            Ub_CEN=Ub_CEN & Regions{TumeurBis}; 
        end
        Dice_CV(i,iter)=dice(Ub_CV,Gold_Standards{i});
        Dice_CP(i,iter)=dice(Ub_CP,Gold_Standards{i});
        Dice_CEN(i,iter)=dice(Ub_CEN,Gold_Standards{i});
        fprintf("fin de l'it�ration %d.%d\n",iter,i);
    end
end

%% Affichage des r�sultats des tests de S�rensen-Dice
f1=figure('Name','Resultats des tests de S�rensen-Dice','NumberTitle','off');
f2=figure('Name','Etude du temps de calcul','NumberTitle','off');
f3=figure('Name','Images','NumberTitle','off');

iter=0;
for i=Im2Train
    iter=iter+1;
    figure(f1);
    
    subplot(1,NbImages,iter);
    Dice=[Dice_CV(i,:);Dice_CEN(i,:);Dice_CP(i,:)];    
    if Log_scale
        semilogx(Lambda,Dice);
    else
        plot(Lambda,Dice);
    end
    legend('CV','CEN','CP','Location','best');
    
    dice_max_CP=max(Dice_CP(i,:));
    tmp1=find(Dice_CP(i,:)==dice_max_CP);
    lambda_max_CP=Lambda(tmp1(1));
    dice_max_CEN=max(Dice_CEN(i,:));
    tmp2=find(Dice_CEN(i,:)==dice_max_CEN);
    lambda_max_CEN=Lambda(tmp2(1));
    dice_max_CV=max(Dice_CV(i,:));
    tmp3=find(Dice_CV(i,:)==dice_max_CV);
    lambda_max_CV=Lambda(tmp3(1));
    
    xlabel({'\lambda'; ['CV : max= ',num2str(dice_max_CV),',  \lambda_{max}= ',num2str(lambda_max_CV)];...
        ['CEN : max= ',num2str(dice_max_CEN),',  \lambda_{max}= ',num2str(lambda_max_CEN)];...
        ['CP : max= ',num2str(dice_max_CP),',  \lambda_{max}= ',num2str(lambda_max_CP)]});
        
    ylabel('Dice index');
%     if i==ceil(NbImages/2)
%         title({"Evolution de l'indice de Dice fonction du param�tre \lambda";['Image ',num2str(i)]});
%     else
%         title(["Image ",num2str(i)]);
%     end
    if i==1
        title('Square');
    elseif i==2
        title('Coins');
    else
        title('Tumor');
    end
    figure(f2)
    subplot(1,NbImages,iter);
    t=reshape(Time(i,:,:),[],length(Lambda));
    if Log_scale
        semilogx(Lambda,t);
    else
        plot(Lambda,t);
    end
    
    legend('CV','CEN','CP','Location','best');
    ylabel('t');
    xlabel('\lambda');
%     if i==ceil(NbImages/2)
%         title({'Evolution du temps de calcul en fonction du param�tre \lambda';['Image ',num2str(i)]});
%     else
%         title(["Image ",num2str(i)]);
%     end
    if i==1
        title('Square');
    elseif i==2
        title('Coins');
    else
        title('Tumor');
    end
    
    figure(f3)
    subplot(1,NbImages,iter)
    imshow(Images{i});
    hold on
    contour(Masks{i},'r','Linewidth',3);
    hold off
    
    title(["Image et masque initial ",num2str(i)]);
end


if ~Import
    save('Images/Masks','Masks');
    save('Images/Backgrounds','Backgrounds');
    save('Images/Regions','Regions');
end