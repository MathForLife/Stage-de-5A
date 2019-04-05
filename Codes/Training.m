clear all; close all;

filename={'Square','GeometricShape','Coins','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Train=[1,2]; 
ImWithRegion=4:7; 

ChangeMasks=false;
Foreground2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3;

[Images, Foregrounds, Backgrounds,Regions,Gold_Standards]=ImportImageMasks(filename,extension,Im2Train,ImWithRegion,ChangeMasks,Foreground2Change,Background2Change,Region2Change);
%% Initialisation des parametres pour le test de Serensen-Dice
fprintf('Initialisation des parametres numeriques\n');
pow_min=-10; pow_max=10; pow_step=1;
Log_scale=true;

pow_length=length(pow_min:pow_step:pow_max);
Lambda=nan(1,pow_length); NbImages=length(Im2Train);

Dice_CV=nan(3,pow_length);
Dice_CP=nan(3,pow_length);
Dice_CEN=nan(3,pow_length);
Time=nan(3,3,pow_length);

%% Parametres numeriques
itermax=500; cinconnu=false;
stop_1=1.e-8; stop_2=1.e-8;
StopConditions=[itermax, stop_1, stop_2];

order=2; eta=1; n=10; epsilon=1;
ParametersCV=[n,order,eta,epsilon];

mu=0.5; tho=1/16;
ParametersCEN=[tho,mu,epsilon];

theta=1; tho_u=0.25; tho_z=0.25;
ParametersCP=[tho_u,tho_z, mu, theta];


fprintf('%d iterations prevues\n',pow_length);
iter=0;

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
    %% Boucle sur les images
    for im=Im2Train
        c1=mean(Images{im}(Foregrounds{im})); c2=mean(Images{im}(~Foregrounds{im}));
        Colors=[c1,c2,cinconnu];
        
        tic
        [Ub_CV, J_CV, err_u,err_J, Niter_CV]=ChanVese(Images{im}, Foregrounds{im}, lambda,  ParametersCV, Colors, StopConditions);
        t1=toc;
        
        tic
        [Ub_CEN, ~, ~,~, ~]=ChanEsedogluNikolova(Images{im}, Foregrounds{im}, lambda, ParametersCEN, Colors, StopConditions);
        t2=toc;
        
        tic
        [Ub_CP, ~, ~,~, ~]=DualFormulation(Images{im}, Foregrounds{im}, lambda, ParametersCP, Colors, StopConditions);
        t3=toc;
        
        Time(im,:,iter)=[t1,t2,t3];
        
        if ismember(im,ImWithRegion)
            Ub_CV=Ub_CV & Regions{im};
            Ub_CP=Ub_CP & Regions{im};
            Ub_CEN=Ub_CEN & Regions{im}; 
        end
        Dice_CV(im,iter)=SorensenDice(Ub_CV,Gold_Standards{im});
        Dice_CP(im,iter)=SorensenDice(Ub_CP,Gold_Standards{im});
        Dice_CEN(im,iter)=SorensenDice(Ub_CEN,Gold_Standards{im});
        fprintf("fin de l'iteration %d.%d\n",iter,im);
    end
end

%% Affichage des resultats des tests de Serensen-Dice
f1=figure('Name','Resultats des tests de Serensen-Dice','NumberTitle','off');
f2=figure('Name','Etude du temps de calcul','NumberTitle','off');
f3=figure('Name','Images','NumberTitle','off');

iter=0;
for im=Im2Train
    iter=iter+1;
    figure(f1);
    
    subplot(1,NbImages,iter);
    Dice=[Dice_CV(im,:);Dice_CEN(im,:);Dice_CP(im,:)];    
    if Log_scale
        semilogx(Lambda,Dice);
    else
        plot(Lambda,Dice);
    end
    legend('CV','CEN','CP','Location','best');
    
    dice_max_CP=max(Dice_CP(im,:));
    tmp1=find(Dice_CP(im,:)==dice_max_CP);
    lambda_max_CP=Lambda(tmp1(1));
    dice_max_CEN=max(Dice_CEN(im,:));
    tmp2=find(Dice_CEN(im,:)==dice_max_CEN);
    lambda_max_CEN=Lambda(tmp2(1));
    dice_max_CV=max(Dice_CV(im,:));
    tmp3=find(Dice_CV(im,:)==dice_max_CV);
    lambda_max_CV=Lambda(tmp3(1));
    
    xlabel({'\lambda'; ['CV : max= ',num2str(dice_max_CV),',  \lambda_{max}= ',num2str(lambda_max_CV)];...
        ['CEN : max= ',num2str(dice_max_CEN),',  \lambda_{max}= ',num2str(lambda_max_CEN)];...
        ['CP : max= ',num2str(dice_max_CP),',  \lambda_{max}= ',num2str(lambda_max_CP)]});
        
    ylabel('Dice coefficient');
    title(filename{im});

    figure(f2)
    subplot(1,NbImages,iter);
    t=reshape(Time(im,:,:),[],length(Lambda));
    if Log_scale
        semilogx(Lambda,t);
    else
        plot(Lambda,t);
    end
    
    legend('CV','CEN','CP','Location','best');
    ylabel('t'); xlabel('\lambda');
    title(filename{im});
    
    figure(f3)
    subplot(1,NbImages,iter)
    imshow(Images{im});
    hold on
    contour(Foregrounds{im},'r','Linewidth',3);
    contour(Backgrounds{im},'g','Linewidth',3);
    hold off
    
    title([filename{im}, " and initial masks"]);
end

%% Sauvegarde des masques modifies
if ChangeMasks
    SaveMasks(filename, Foregrounds,Foreground2Change,Backgrounds, Background2Change,Regions,Region2Change)
end