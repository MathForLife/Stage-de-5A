clear all; close all;

addpath(genpath('../Images/'))
load('Gold_Standards.mat');
load('Images/Masks','Masks');
%load('Images/Backgrounds','Backgrounds');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Affiche des masques
% coins_mask=Gold_Standards{2};
% phi=signed_distance_from_mask(coins_mask);
% surface(phi,'EdgeColor','none');
% temp=abs(phi)==1;
% hold on;
% contour3(phi.*temp,'r');
% hold off;
% colorbar();
% figure;
% imshow(double(coins_mask));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Affiche la fonction de Heaviside
% n=2;
% x=linspace(-5,5,100);
% H=zeros(size(x)); H(x>0)=1;
% y1=Heavyside_eta(x,1,2);
% y2=Heavyside_eta(x,0.5,2);
% y3=Heavyside_eta(x,0.1,2);
% Y=[y3;y2;y1];
% figure;
% plot(x,H,'-b');
% ylim([-0.5,1.5]);
% xlabel("x"); ylabel("H(x)");
% figure;
% plot(x,Y)
% legend('\eta = 0.1','\eta = 0.5','\eta = 1');
% xlabel("x"); ylabel("H_{\eta}(x)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Affichage des résultats pour les différentes méthodes
% lambda_CP=100; lambda_CEN=100; lambda_CV=1000;
% Import=true;
% Bruitage=false;
% Im2Test=[3]; Mask2Change=2; Background2Change=2;
% % Parametres numeriques
% itermax=500; cinconnu=false;
% mu=0.5; eps=1; theta=1;
% tho_u=0.5; tho_z=0.25;
% stop_1=-1.e-6; stop_2=-1.e-6;
% order=2; n=10; eta=1; %Paramètres nécessaires à CV
% % Lecture des images
% Image1=double(imread('Square.jpg'));
% Image2=double(imread('eight.tif'));
% Image3=double(imread('TumeurCerveaubis.png'));
%
% NbImages=length(Im2Test);
% Images={Image1,Image2,Image3};
% % Normalisation des images
% for i=Im2Test
%     Images{i}=Image_Normalisation(Images{i},"2D");
% end
% %% Ajout de bruit aux images
% if Bruitage==true
%     bruit=0.1;
%
%     for i=Im2Test
%         Images{i}=Images{i}+bruit*randn(size(Images{i}));
%         Images{i}=Image_Normalisation(Images{i},"2D");
%     end
% end
%
% % Importation des masques ou creation de ces derniers
% if ~Import
%     for m=Mask2Change
%         Masks{m}=roipoly(Images{m});
%     end
%     close;
%     for b=Background2Change
%         Backgrounds{b}=roipoly(Images{b});
%     end
%     close;
% end
%
% % Creation des figures
% for i=Im2Test
%     c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(Backgrounds{i}));
%     [Ub_CV, J_CV, err_u,err_J, Niter_CV]=ChanVese(Images{i}, double(Masks{i}), lambda_CV, eps, eta,stop_1,stop_2,c1,c2, cinconnu, n, order, itermax);
%     [Ub_CP, J_CP, Err_u_Cp,Err_J_CP, Niter_CP]=DualFormulation(Images{i}, double(Masks{i}), lambda_CP, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
%     [Ub_CEN, J_CEN, Err_u_CEN,Err_J_CEN, Niter_CEN]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda_CEN, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
%
%     figure;
%     subplot(1,3,1);
%     imagesc(Images{i}); axis off; axis image;
%     colormap gray
%     hold on
%     contour(Ub_CV,'r','Linewidth',3);
%     hold off
%     title({'Chan-Vese';['\lambda_{CV}= ',num2str(lambda_CV)]})
%     xlabel(['\lambda_CV= ',num2str(lambda_CV)]);
%
%     subplot(1,3,2);
%     imagesc(Images{i}); axis off; axis image;
%     colormap gray
%     hold on
%     contour(Ub_CEN,'r','Linewidth',3);
%     hold off
%     title({'Chan-Esedoglu-Nikolova';['\lambda_{CEN}= ',num2str(lambda_CEN)]})
%     xlabel(['\lambda_CEN= ',num2str(lambda_CEN)]);
%
%     subplot(1,3,3);
%     imagesc(Images{i}); axis off; axis image;
%     colormap gray
%     hold on
%     contour(Ub_CP,'r','Linewidth',3);
%     hold off
%     title({'Chambolle-Pock';['\lambda_{CP}= ',num2str(lambda_CP)]})
%
% %
% %     figure;
% %     imagesc(Images{i}); axis off; axis image;
% %     colormap gray
% %     hold on
% %     contour(Masks{i},'r','Linewidth',3);
% %     contour(Backgrounds{i},'g','Linewidth',3);
% %     hold off
% %     legend('Initial segmentation mask','Background mask');
% %     legend('Initial Foreground','Initial Background');
%
%     figure;
%     subplot(3,1,1)
%     plot(1:Niter_CV,J_CV); xlabel('iterations'); ylabel('J(u_b)');
%     title('CV');
%     subplot(3,1,2)
%     plot(1:Niter_CEN,J_CEN); xlabel('iterations'); ylabel('J(u_b)');
%     title('CEN');
%     subplot(3,1,3)
%     plot(1:Niter_CP,J_CP); xlabel('iterations'); ylabel('J(u_b)');
%     title('CP');
%
% end
%% Normalisation des differents masques
filename={'Square','GeometricShape','Coins','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
Im2Def=3;
Foregrounds={};
Backgrounds={};
Regions={}; %a=figure; b=figure; c=figure;
for i=Im2Def
%     load(['Foregrounds/',filename{i},'_FG.mat'],'FG');
%     Foregrounds{i}=FG;
%     load(['Backgrounds/',filename{i},'_BG.mat'],'BG');
%     Backgrounds{i}=BG;
%     figure(a)
%     imagesc(FG);
%     figure(b)
%     imagesc(BG);
%     if i>3
%         load(['../Images/Regions/',filename{i},'_Region.mat'],'R');
%         Regions{i}=R;
%         figure(c);
%         imagesc(R);
%     end
    FG=double(imread(['Foregrounds/',filename{i},'_FG.png']));
    FG=Image_Normalisation(FG,'2D');
    FG(FG<1)=0; FG=logical(FG);
    save(['../Images/Foregrounds/',filename{i},'_FG.mat'],'FG');
    
    BG=double(imread(['Backgrounds/',filename{i},'_BG.png']));
    BG=Image_Normalisation(BG,'2D');
    BG(BG<1)=0; BG=logical(BG);
    save(['../Images/Backgrounds/',filename{i},'_BG.mat'],'BG');
    
    GS=double(imread(['Gold_Standards/',filename{i},'_GS.png']));
    GS=Image_Normalisation(GS,'2D');
    GS(GS<1)=0; GS=logical(GS);
    save(['../Images/Gold_Standards/',filename{i},'_GS.mat'],'GS');
    if strcmp(filename{i},'Coins')
        FG=double(imread(['Foregrounds/',filename{i},'_FG_Glob.png']));
        FG=Image_Normalisation(FG,'2D');
        FG(FG<1)=0; FG=logical(FG);
        save(['../Images/Foregrounds/',filename{i},'_FG_Glob.mat'],'FG');
        
        BG=double(imread(['Backgrounds/',filename{i},'_BG_Glob.png']));
        BG=Image_Normalisation(BG,'2D');
        BG(BG<1)=0; BG=logical(BG);
        save(['../Images/Backgrounds/',filename{i},'_BG_Glob.mat'],'BG');
    end
    
    if i>3
        R=double(imread(['../Images/Regions/',filename{i},'_Region.png']));
        R=Image_Normalisation(R,'2D');
        R(R<1)=0; R=logical(R);
        save(['../Images/Regions/',filename{i},'_Region.mat'],'R');
    end
end