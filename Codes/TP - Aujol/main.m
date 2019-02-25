%% Initialisation 
addpath(genpath('../../'));
Image=double(imread('TumeurCerveaubis.png'));
Image=Image(:,:,1);

bruitage=false;
if bruitage==true
    bruit=0.05;
    Image=Image+bruit*255*randn(size(Image));
    
    Min=min(min(min(Image)),0);
    Max=max(max(max(Image)),255);
    
    Image=Image-Min;
    Image=Image*255/(Max-Min);
end
mask=roipoly(Image/255.);
imshow(mask);
% Definition des parametres numeriques 
lambda=1.e-4;
eta=2;
eps=1;
sigma=0.1;

reset_LS=5; % Nombre d'iterations de descente de gradient avant de recalculer la level-set
choice_heavyside=2; %1 : regularisation C^2  , 2 : regularisation C^inf (cf. Active Contours Whithout Edges) 
% Remarque : si eta<2, le choix 1 va creer une mesure de Dirac nulle partout, donc inutilisable.
itermax=200;

%% Algorithme de Chan - Vese
tic
[phib,JCV, TimeStepCV, StopCondCV, c1, c2, niterCV]=chanvese(Image, mask, lambda, sigma, eps, eta, reset_LS, choice_heavyside, itermax);
t1=toc;

% figure();
% subplot(311);
% %semilogy(1:niter,CostFunc)
% plot(1:niterCV,JCV)
% title('Evolution de la fonction cout')
% xlabel('iterations')
% subplot(312);
% plot(1:niterCV,TimeStepCV);
% title('Evolution du pas de descente du gradient')
% xlabel('iterations')
% subplot(313);
% plot(1:niterCV,StopCondCV);
% title("Evolution du critere d'arret")
% xlabel('iterations')

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(phib,'r','Linewidth',3);
title(['Resultat Chan-Vese pour ' num2str(niterCV),' iterations']);
sprintf('C1 = %5.1f, C2= %5.1f, niter= %d',c1,c2,niterCV)


%% Algorithme de Chan-Esedoglu-Nikolova
mu=0.5;
tho=0.1;

c1=110;
c2=227;

cinconnu=true;

tic
[ubCEN, JCEN, StopCondCEN, niterCEN]=ChanEsedogluNikolova(Image, mask, lambda, mu, tho, sigma, eps, c1, c2, cinconnu, itermax);
t2=toc;

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(ubCEN,'r','Linewidth',3);
title(['Resultat Chan-Esedoglu-Nikolova pour ',num2str(niterCEN),' iterations']);

%% Algorithme Primal-Dual
lambda=1.e-4;
tho=0.25;
theta=0;
Projstep=0.5;

tic
[ubPD, J1PD, J2PD, StopCondPD, niterPD]=DualFormulation(Image, mask, lambda, mu, tho, theta, Projstep, sigma, c1, c2, cinconnu, itermax);
t3=toc;

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(ubPD,'r','Linewidth',3);
title(['Resultat Primal-Dual pour ',num2str(niterPD),' iterations']);

%% Algorithme Primal-Dual acceléré 
theta=1;

tic
[ubPDA, J1PDA, J2PDA, StopCondPDA, niterPDA]=DualFormulation(Image, mask, lambda, mu, tho, theta, Projstep, sigma, c1, c2, cinconnu, itermax);
t4=toc;

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(ubPD,'r','Linewidth',3);
title(['Resultat Accelerated-Primal-Dual pour ',num2str(niterPDA),' iterations']);

%% Courbes de convergence
figure();
subplot(221);
%semilogy(1:niter,CostFunc)
plot(1:niterCV,JCV)
title('Fonction cout algo Chan-Vese')
xlabel('iterations')
subplot(222);
%semilogy(1:niterCEN,J)
plot(1:niterCEN,JCEN);
title('Fonction cout algo Chan-Esedoglu-Nikolova')
xlabel('iterations')
subplot(223);
%semilogy(1:niterCEN,J)
plot(1:niterPD,J1PD);
title('Fonction cout algo Primal-Dual')
xlabel('iterations')
subplot(224);
%semilogy(1:niterCEN,J)
plot(1:niterPDA,J1PDA);
title('Fonction cout algo Accelerated-Primal-Dual')
xlabel('iterations')

sprintf("Temps d'execution des algorithmes :\n Chan-Vese : %6.4f",t1)
sprintf('Chan-Esedoglu-Nikolova : %6.4f ',t2)
sprintf('Primal-Dual : %6.4f ',t3)
sprintf('Accelerated-Primal-Dual : %6.4f ',t4)