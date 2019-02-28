addpath(genpath('../../'));
Image1=zeros(248);
Image1(62:186,62:186)=255;
Image2=double(imread('eight.tif'));
Image3=double(imread('TumeurCerveaubis.png'));
Image3=Image3(:,:,1);

%% Initialisation des booléans
Images={Image1,Image2,Image3};
Plot=false;
Import=true;
Dice_test=true;
%% Normalisation des images
for i=1:3
    Images{i}=Image_Normalisation(Images{i},"2D");
end
%% Ajout de bruit aux images
bruitage=false;
if bruitage==true
    bruit=0.05;
    
    for i=1:3
        Images{i}=Images{i}+bruit*255*randn(size(Images{i}));
        Images{i}=Image_Normalisation(Images{i},"2D");
    end
end

if Import
    load('Images/Masks','Masks');
else
    mask1=roipoly(Images{1});
    mask2=roipoly(Images{2});
    mask3=roipoly(Images{3}); close;
    Masks={mask1,mask2,mask3};
end
%% Paramètres numériques
itermax=200; cinconnu=false;
mu=0.1; eps=1;
theta=1;
tho_u=0.5; tho_z=0.25;
stop_1=1.e-6; stop_2=1.e-6;

%% Initialisation des paramètres pour le test de Sørensen-Dice 
pow_min=1; pow_max=2;
if Dice_test
    load('Gold_Standards.mat');
    
    Lambda=nan(1,abs(pow_max-pow_min+1));
    Dice_CP=nan(3,abs(pow_max-pow_min+1));
    Dice_CEN=nan(3,abs(pow_max-pow_min+1));
    fprintf('%d itérations prévues\n',length(Lambda));
    iter=0;
end
%% Boucle des lambdas
%for ind=pow_min:pow_max
for ind=1:10:500
    %lambda=ind*10^1;
    lambda=ind;
    if Plot
        f1=figure('Name','Resultat de la segmentation','NumberTitle','off'); 
        f2=figure('Name','Evolution des quantités de contrôle','NumberTitle','off');
    elseif Dice_test
        iter=iter+1;
        fprintf("Début de l'itération %d\n",iter);
    end
    
    for i=1:3
        c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(~Masks{i}));
        
        if Plot
            [Ub, J, Err_u,Err_J, Niter]=DualFormulation(Images{i}, double(Masks{i}), lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
            %[Ub, J, Err_u,Err_J, Niter]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
        
            figure(f1);
            subplot(1,3,i)
            imagesc(Images{i}); axis off; axis image;
            colormap gray
            hold on
            contour(Ub,'r','Linewidth',3);
            title(['Segmentation par C-P avec \lambda= ',num2str(lambda),' et pour ',num2str(Niter),' iterations']);
            hold off

            figure(f2);
            subplot(3,3,i);
            plot(1:Niter,J)
            title(['Fonction cout pour \lambda= ',num2str(lambda)])
            xlabel('iterations')
            subplot(3,3,i+3)
            semilogy(1:Niter,Err_u)
            title('Erreur relative sur u_{n+1}-u_n')
            xlabel('iterations')
            subplot(3,3,i+6)
            semilogy(1:Niter,Err_J)
            title('Erreur relative sur J(u_{n+1})-J(u_n)')
            xlabel('iterations')
        elseif Dice_test
            Lambda(iter)=lambda;
            [Ub_CP, ~, ~,~, ~]=DualFormulation(Images{i}, double(Masks{i}), lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
            [Ub_CEN, ~, ~,~, ~]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
            
            if i==3
                Dice_CP(i,iter)=dice(Ub_CP&mask_quantif,Gold_Standards{i}&mask_quantif);
                Dice_CEN(i,iter)=dice(Ub_CEN&mask_quantif,Gold_Standards{i}&mask_quantif);
                fprintf("fin de l'itération %d.%d\n",iter,i);
            else
                Dice_CP(i,iter)=dice(Ub_CP,Gold_Standards{i});
                Dice_CEN(i,iter)=dice(Ub_CEN,Gold_Standards{i});
                fprintf("fin de l'itération %d.%d\n",iter,i);
            end
        end
    end
end

%% Affichage des résultats des tests de Sørensen-Dice
if Dice_test && ~Plot
    figure('Name','Resultats des tests de Sørensen-Dice','NumberTitle','off');
    for i=1:3
        subplot(2,3,i);
        %semilogx(Lambda,Dice_CP(i,:));
        plot(Lambda,Dice_CP(i,:));
        xlabel('\lambda');
        ylabel('Dice index');
        title(["Chambol-Pock sur l'image ",num2str(i)]);
        subplot(2,3,i+3);
        %semilogx(Lambda,Dice_CEN(i,:));
        plot(Lambda,Dice_CEN(i,:));
        xlabel('\lambda');
        ylabel('Dice index');
        title("Chan-Esedoglu-Nikolova");
        fprintf('lambda_{max} : %4.2f, value : %4.2f',Lambda(Dice_CP(3,:)==max(Dice_CP(3,:))),max(Dice_CP(3,:)));
    end
end

if ~Import
    save('Images/Masks','Masks');
end