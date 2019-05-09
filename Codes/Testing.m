clear all; close all;

image_names={'Square','GeometricShape','Coins','Flag','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Test=[5]; NbImages=length(Im2Test);
ImWithRegion=4:7; 

ChangeMasks=false; Bruitage=false; cinconnu=false;   ChooseTexture=true;
Foreground2Change=[5]; Background2Change=[]; Region2Change=[];

[Images, Textures, Foregrounds, Backgrounds,~,~]=ImportImageMasks(image_names,extension,Im2Test,ImWithRegion,ChooseTexture,ChangeMasks,Foreground2Change,Background2Change,Region2Change);

%% Parametres numeriques
fprintf('Initialisation des parametres \n');
lambda_CV=1.e2; lambda_CEN=1.e3; lambda_CP=1.e2;

itermax=500; 
stop_1=-1.e-8; stop_2=-1.e-8;
StopConditions=[itermax, stop_1, stop_2];

order=2; eta=1; n=10; epsilon=1;
ParametersCV=[n,order,eta,epsilon];

mu=0.5; tho=1/sqrt(8);
ParametersCEN=[tho,mu,epsilon];

theta=1; tho_u=0.25; tho_z=0.25;
ParametersCP=[tho_u,tho_z, mu, theta];

%% Ajout de bruit aux images
if Bruitage==true
    bruit=0.05;
    
    for im=Im2Test
        Images{im}=Images{im}+bruit*randn(size(Images{im}));
        Images{im}=Image_Normalisation(Images{im},"2D");
    end
end

%% Definition des figures
f1=figure('Name','Evolution des quantites de controle pour C-V','NumberTitle','off');
f2=figure('Name','Evolution des quantites de controle pour C-E-N','NumberTitle','off');
f3=figure('Name','Evolution des quantites de controle pour C-P','NumberTitle','off');
f4=figure('Name','Images','NumberTitle','off');

%% Boucle sur les images
iter=0;
for im=Im2Test
    iter=iter+1;
    
%     c1=mean(Images{im}(Foregrounds{im})); c2=mean(Images{im}(Backgrounds{im}));
%     Colors=[c1,c2, cinconnu];
    
%     [Ub_CV, J_CV, Err_u_CV,Err_J_CV, Niter_CV]=ChanVese(Images{im}, Foregrounds{im}, lambda_CV, ParametersCV, Colors, StopConditions);
%     [Ub_CEN, J_CEN, Err_u_CEN,Err_J_CEN, Niter_CEN]=ChanEsedogluNikolova(Images{im}, Foregrounds{im}, lambda_CEN, ParametersCEN, Colors, StopConditions);
%     [Ub_CP, J_CP, Err_u_Cp,Err_J_CP, Niter_CP]=DualFormulation(Images{im}, Foregrounds{im}, lambda_CP, ParametersCP,Colors, StopConditions);
%     
    [Ub_CV, J_CV, Err_u_CV,Err_J_CV, Niter_CV]=ChanVese_Texture(Images{im}, Textures{im}, lambda_CV, Foregrounds{im}, Backgrounds{im}, cinconnu, ParametersCV, StopConditions);
    [Ub_CEN, J_CEN, Err_u_CEN,Err_J_CEN, Niter_CEN]=ChanEsedogluNikolova_Texture(Images{im}, Textures{im}, lambda_CEN, Foregrounds{im}, Backgrounds{im}, cinconnu, ParametersCEN, StopConditions);
    [Ub_CP, J_CP, Err_u_Cp,Err_J_CP, Niter_CP]=DualFormulation_Texture(Images{im}, Textures{im}, lambda_CP, Foregrounds{im}, Backgrounds{im}, cinconnu, ParametersCP, StopConditions);
    
    figure;
    subplot(1,3,1)
    imagesc(Images{im}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CV,'r','Linewidth',3);
    hold off
    title({'Chan-Vese';['\lambda_{CV}= ',num2str(lambda_CV)]});
    
    subplot(1,3,2)
    imagesc(Images{im}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CEN,'r','Linewidth',3);
    hold off
    title({'Chan-Esedoglu-Nikolova';['\lambda_{CEN}= ',num2str(lambda_CEN)]});
    
    subplot(1,3,3)
    imagesc(Images{im}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CP,'r','Linewidth',3);
    hold off
    title({'Chambolle-Pock';['\lambda_{CP}= ',num2str(lambda_CP)]});
    
    figure(f1);
    subplot(3,NbImages,iter);
    plot(J_CV)
    title({image_names{im};['Fonction cout pour \lambda= ',num2str(lambda_CV)]})
    xlabel('iterations')
    subplot(3,NbImages,iter+NbImages)
    semilogy(Err_u_CV)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,NbImages,iter+2*NbImages)
    semilogy(Err_J_CV)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')
    
    figure(f2);
    subplot(3,NbImages,iter);
    plot(J_CEN)
    title({image_names{im};['Fonction cout pour \lambda= ',num2str(lambda_CEN)]})
    xlabel('iterations')
    subplot(3,NbImages,iter+NbImages)
    semilogy(Err_u_CEN)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,NbImages,iter+2*NbImages)
    semilogy(Err_J_CEN)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')
    
    figure(f3);
    subplot(3,NbImages,iter);
    plot(J_CP)
    title({image_names{im};['Fonction cout pour \lambda= ',num2str(lambda_CP)]})
    xlabel('iterations')
    subplot(3,NbImages,iter+NbImages)
    semilogy(Err_u_Cp)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,NbImages,iter+2*NbImages)
    semilogy(Err_J_CP)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')    
    
    figure(f4)
    subplot(1,NbImages,iter)
    imshow(Images{im});
    hold on
    contour(Foregrounds{im},'r','Linewidth',3);
    contour(Backgrounds{im},'g','Linewidth',3);
    hold off
    legend('\partial\Omega_1','\partial\Omega_0')
    title([image_names{im}, " and initial masks"]);
end
