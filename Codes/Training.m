addpath('../')

%% Importation des gold standards + var booléennes
load('Gold_Standards.mat');
Import=true;
Log_scale=true;

%% Lecture des images
Image1=zeros(248);
Image1(62:186,62:186)=255;
Image2=double(imread('eight.tif'));
Image3=double(imread('TumeurCerveaubis.png'));

NbImages=3;  %Nombre d'images testées
TumeurBis=3; %Numéro de l'image contenant le détail de la tumeure
Images={Image1,Image2,Image3};
%% Normalisation des images
for i=1:NbImages
    Images{i}=Image_Normalisation(Images{i},"2D");
end

% Importation des masques ou création de ces derniers
if Import
    load('Images/Masks','Masks');
    load('Images/Region_tumeur_bis','Region_tumeur_bis')
else
    mask1=roipoly(Images{1});
    mask2=roipoly(Images{2});
    mask3=roipoly(Images{3}); close;
    Masks={mask1,mask2,mask3};
    
    Region_tumeur_bis=roipoly(Images{3}); close;
end

%% Paramètres numériques
itermax=200; cinconnu=false;
mu=0.1; eps=1;
theta=1;
tho_u=0.5; tho_z=0.25;
stop_1=1.e-6; stop_2=1.e-6;

%% Initialisation des paramètres pour le test de Sørensen-Dice
pow_min=0; pow_max=5; pow_step=1;

Lambda=nan(1,abs(pow_max-pow_min+1));
Time=nan(3,2,abs(pow_max-pow_min+1));
Dice_CP=nan(3,abs(pow_max-pow_min+1));
Dice_CEN=nan(3,abs(pow_max-pow_min+1));
fprintf('%d itérations prévues\n',length(Lambda));
iter=0;

%% Boucle des lambdas
for ind=pow_min:pow_step:pow_max
    if Log_scale
        lambda=10^ind;
    else
        lambda=ind;
    end
    iter=iter+1;
    fprintf("Début de l'itération %d\n",iter);
    
    Lambda(iter)=lambda;
    
    for i=1:NbImages
        c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(~Masks{i}));
        tic
        [Ub_CP, ~, ~,~, ~]=DualFormulation(Images{i}, double(Masks{i}), lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
        t1=toc;
        
        tic
        [Ub_CEN, ~, ~,~, ~]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
        t2=toc;
        
        if i==TumeurBis
            Ub_CP=Ub_CP & Region_tumeur_bis;
            Ub_CEN=Ub_CEN & Region_tumeur_bis;
        end
        Dice_CP(i,iter)=dice(Ub_CP,Gold_Standards{i});
        Dice_CEN(i,iter)=dice(Ub_CEN,Gold_Standards{i});
        Time(i,:,iter)=[t1,t2];
        fprintf("fin de l'itération %d.%d\n",iter,i);
    end
end


%% Affichage des résultats des tests de Sørensen-Dice
f1=figure('Name','Resultats des tests de Sørensen-Dice','NumberTitle','off');
f2=figure('Name','Etude du temps de calcul','NumberTitle','off');
f3=figure('Name','Images','NumberTitle','off');
for i=1:NbImages
    figure(f1);
    subplot(1,3,i);
    Dice=[Dice_CP(i,:);Dice_CEN(i,:)];
    if Log_scale
        semilogx(Lambda,Dice);
    else
        plot(Lambda,Dice);
    end
    legend('CP','CEN','Location','best');
    
    dice_max_CP=max(Dice_CP(i,:));
    tmp1=find(Dice_CP(i,:)==dice_max_CP);
    lambda_max_CP=Lambda(tmp1(1));
    dice_max_CEN=max(Dice_CEN(i,:));
    tmp2=find(Dice_CEN(i,:)==dice_max_CEN);
    lambda_max_CEN=Lambda(tmp2(1));
    
    xlabel({'\lambda'; ['CP : max= ',num2str(dice_max_CP),',  \lambda_{max}= ',num2str(lambda_max_CP)];...
        ['CEN : max= ',num2str(dice_max_CEN),',  \lambda_{max}= ',num2str(lambda_max_CEN)]});
    ylabel('Dice index');
    if i==ceil(NbImages/2)
        title({"Evolution de l'indice de Dice fonction du paramètre \lambda";['Image ',num2str(i)]});
    else
        title(["Image ",num2str(i)]);
    end
    
    figure(f2)
    subplot(1,3,i);
    t=reshape(Time(i,:,:),2,length(Lambda));
    if Log_scale
        semilogx(Lambda,t);
    else
        plot(Lambda,t);
    end
    legend('CP','CEN','Location','best');
    ylabel('t');
    xlabel('\lambda');
    if i==ceil(NbImages/2)
        title({'Evolution du temps de calcul en fonction du paramètre \lambda';['Image ',num2str(i)]});
    else
        title(["Image ",num2str(i)]);
    end
    
    figure(f3)
    subplot(1,3,i)
    imshow(Images{i});
    title(["Image ",num2str(i)]);
end


if ~Import
    save('Images/Masks','Masks');
    save('Images/Region_tumeur_bis','Region_tumeur_bis')
end