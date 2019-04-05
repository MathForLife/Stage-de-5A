function [phib, J, err_u,err_J, niter]=ChanVese(Image, mask, lambda, Parameters, Colors, StopConditions)
%% Algorithme de Chan-Vese calculant un masque binaire de segmentation a partir d'une image et d'un masque initial
% INPUTS:
% Image: Image que l'on cherche à segmenter en 2 régions. L'image est en nuances de gris et peut etre en 2D ou en 3D
% mask : Masque représentant une partie de la zone que l'on souhaite segmenter. Il est possible de faire une union de masques
% lambda: Parametre de controle entre attache aux donnees et regularisation sur la fonction J à minimiser
% Parameters: Contient des informations propres a l'implementation numerique de CV
%     -1 n : nombre d'iterations avant de reinitialiser la fonction Level-Set
%     -2 order : indique quelle methode de lissage est utilisee pour les
%     fonctions de Heaviside et de Dirac (peut prendre la valeur 1 ou 2)
%     -3 eta : parametre de lissage pour les fonctions de Heaviside et de Dirac 
%     -4 epsilon: Parametre de lissage pour le calcul de la norme 
% Colors: Contient des informations sur les parametres de couleurs donnes utilises pour la minimisation
%     -1 c1 : couleur moyenne de la region a segmenter
%     -2 c2 : couleur moyenne du fond de l'image
%     -3 cinconnu : variable booleenne indiquant si on modifie les couleurs a chaque iteration (influe sur le caractere convexe de la fonction J)
% StopConditions: Contient les conditions d'arrêt de l'algorithme
%     -1 itermax: Maximum d'itérations dans la boucle while
%     -2 stop_u: Condition sur la décroissance de la fonction de masquage u
%     -3 stop_J: Condition sur la décroissance de la fonctionnelle J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% phib: Masque binaire sur l'image de la région segmentée
% J : Vecteur contenant l'evolution de la fonction cout au cours des itérations
% err_u: Vecteur contenant l'evolution du premier critère d'arrêt au cours des itérations
% err_J: Vecteur contenant l'evolution du deuxième critère d'arrêt au cours des itérations
% niter: Nombre d'itérations total de l'algorithme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut algo
itermax=StopConditions(1); c1=Colors(1); c2=Colors(2);
n=Parameters(1); order=Parameters(2); eta=Parameters(3); epsilon=Parameters(4);
%% Initialisation de toutes les variables
niter=1;

J=nan(1,itermax);
err_u=nan(1,itermax);
err_J=nan(1,itermax);

% Initialisation de la fonction Level-Set
phi=signed_distance_from_mask(mask);

% Calcul des fonctions de Heaviside et de Dirac lissees
Hn=Heavyside_eta(phi, eta, order);
dn=delta_eta(phi,eta, order);

I1=(Image-c1).^2;
I2=(Image-c2).^2;

% Calcul des quantités de controle
J(1)=compute_energy_smooth(mask,I1,I2,lambda,epsilon);
cond_u=10; err_u(1)=cond_u; 
cond_J=10; err_J(1)=cond_J;

%% Boucle while
while (niter<itermax && cond_u>StopConditions(2)&& cond_J>StopConditions(3))   
    niter=niter+1;
    
    Hold=Hn; dold=dn;
    
    if Colors(3)
        % Cas ou les couleurs sont aussi des variables de la fonctionnelle
        c1=sum(sum(Image.*Hold))/sum(sum(Hold));
        c2=sum(sum(Image.*(1-Hold)))/sum(sum(1-Hold));
        
        I1=(Image-c1).^2;
        I2=(Image-c2).^2;
    end
    
    % Calcul du gradient de phi et de sa norme lissee
    phix=gradx(phi);
    phiy=grady(phi);
    n_eps=norm_eps(phix,phiy,epsilon);
    
    GradJ=dold.*(div(phix./n_eps,phiy./n_eps)-lambda*(I1-I2));
    
    tho=0.5/max(max(GradJ)); % Ce choix de pas adaptatif permet d'optimiser la descente de gradient
    phi=phi+tho*GradJ;
    phib=phi>0;
    
    if mod(niter,n)==0
        % Au bout de n iterations, on redefinie la fonction Level-Set comme
        % une fonction de distance signee au masque phib
        phi=signed_distance_from_mask(phib);
    end
    
    Hn=Heavyside_eta(phi, eta, order);
    dn=delta_eta(phi,eta, order);
    
    % Calcul des quantitees de controle
    J(niter)=compute_energy_smooth(phib,I1,I2,lambda,epsilon);

    cond_u=norm(Hn-Hold)/norm(Hold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    err_u(niter)=cond_u;
    err_J(niter)=cond_J;
end
%% Redimentionnement des vecteurs de controle de l'algo
err_u(isnan(err_u))=[];
err_J(isnan(err_J))=[];
J(isnan(J))=[];
end