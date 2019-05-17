%% Image segmentation using histograms (if using please refer to the works below).
%% Require IP toolbox
%% By Nicolas Papadakis
%% Multi-class implementation of
%% R.Yildizoglu, J-F. Aujol and N. Papadakis, A convex formulation for global histogram based binary segmentation",  EMMCVPR 2013
%% following the scheme of
%% N. Papadakis and J. Rabin, Convex Histogram-Based Joint Image Segmentation with Regularized Optimal Transport Cost, JMIV, 2017
%% This is a l1 bin to bin metric between color histograms, the Wasserstein version is not implemented here.
%% This example will segment the image im/zebra.png in nb_label=2 regions, 
%% using as reference histograms the color histograms of the regions given in the 2 image masks: im/zebra_scrib1.png and im/zebra_scrib2.png.
%% The histograms of each color channel is discretized in nb_bin bins. The regularization parameter is rho.

close all
clear all
clc
addpath('tools');
%parameters:
list_im = {'zebra', 'aras','chat','BrainTumor','BrainTumorDetail','BrainMeta_A','GBM_A','Gliome003_S','Parenchyme_C','pneumopath_6_A'};  %image to_load
im_file = list_im{5}; %
type='.png';%'.png'; %or '.jpg'


rho=5;   %TV regularization (6 for parrots, 2 for zebra) 
nb_bin=10;  %number of histogram bin on each canal
nb_label=2; %number of class  (given input: 3 for parrots, 2 for zebra)
manual_input_histogram_areas=0; %set 0 if reading input areas
gradient_histograms=false;  %set 1 to consider histograms of gradient norms instead of histograms of colors
niter_max=500.;
gray_image=false;   %set to 1 to force grayvalues
nb_label=max(2,nb_label);
disp(['Segmentation of the image in ' num2str(nb_label) ' regions.']);





if gradient_histograms
    disp('Segment w.r.t histogram of gradient norms on each color channel');
else
    disp('Segment w.r.t histogram of colors');
end


%% load image and masks
disp(['Load image']);
% read image
I=double(imread(['im/' im_file ,type]));
nx=size(I,1);
ny=size(I,2);

%Creation of reference histograms form image masks

Mask_ref=zeros(nx,ny,nb_label);
if manual_input_histogram_areas %select manulally the region to compute the reference histogram
        disp('Choose manually the regions were to compute the reference histograms');
    h=figure;
    for k=1:nb_label
        figure(h)
        imagesc(I/255);
        text_title=['Select the area for the reference histogram of class ' num2str(k)];
        
        set(gcf,'Name',text_title,'NumberTitle','off');
        
        Mask_ref(:,:,k)=roipoly();        
    end
    close(h)
    drawnow;
else %read the region from a file
    disp(['Load masks, one file is required for each of the ' num2str(nb_label) ' regions']);
    for k=1:nb_label        
        M=imread(['im/' im_file '_scrib' num2str(k) ,type]);
        Mask_ref(:,:,k) = M(:,:,1);
    end
end

nx=size(I,1);
ny=size(I,2);
nz=size(I,3);

if gray_image
    nz=1;
end



%% create color histograms
disp(['Creation of input histograms and costs']);

nb_pixel=nx*ny;

nb_bin_total=nb_bin^nz;
histo_ref=sparse(zeros(nb_bin_total,nb_label));

histoIm=sparse(zeros(nb_bin_total,1));

Ifeat=I;
max_feat=256*ones(nz,1);
if gradient_histograms
    Ifeat=I*0;
    for k=1:nz   
        Ifeat(:,:,k)=sqrt(gradx(I(:,:,k)).^2+grady(I(:,:,k)).^2);
        max_feat(k)=max(max(Ifeat(:,:,k)))+1.;
    end
end
    
    
    

for i=1:nx
    for  j=1:ny  
        index=1;
        for k=1:nz
            c=uint16(floor((Ifeat(i,j,k)/max_feat(k))*nb_bin))+1;
            index=index+(c-1)*nb_bin^(nz-k);
        end
%         c1=uint16(floor((Ifeat(i,j,1)/max_feat(1))*nb_bin))+1;
%         c2=uint16(floor((Ifeat(i,j,2)/max_feat(2))*nb_bin))+1;
%         c3=uint16(floor((Ifeat(i,j,3)/max_feat(3))*nb_bin))+1;
%         
%         index=(c1-1)*nb_bin^2+(c2-1)*nb_bin+c3
        
        histoIm(index)=histoIm(index)+1;
        for k=1:nb_label
            if Mask_ref(i,j,k)>0
                histo_ref(index,k)=histo_ref(index,k)+1;               
            end        
        end       
    end
end


histo_ref=(double(histo_ref));
histoIm=(double(histoIm));

indexIm=(histoIm(:)>0);
nb_bin_totalIm=sum(indexIm(:))+1-1;

for k=1:nb_label
    index_ref{k}=(histo_ref(:,k)>0);
end


%bin grid must be the same for bin to bin comparisons
for k=1:nb_label
    index_ref{k}= indexIm;
end

nb_bin_total_ref=zeros(nb_label,1);
for k=1:nb_label   
    nb_bin_total_ref(k)=sum(index_ref{k})+1-1;%to desparsify
    H_ref{k}=zeros(nb_bin_total_ref(k),1);
    H_ref{k}=double(histo_ref(index_ref{k},k));
    H_ref{k}=double(H_ref{k})/sum(H_ref{k}(:));

    %Mean color of bin to create H and cost matrices
    color_ref{k}=zeros(nb_bin_total_ref(k),nz);
end

colorIm=zeros(nb_bin_totalIm,nz);
cptIm=1;
cpt_ref=ones(nb_label,1);

for i=1:nb_bin^nz
    for k=1:nb_label
        if index_ref{k}(i)==1           
            cpt_ref(k)=cpt_ref(k)+1;           
        end
    end
    if indexIm(i)==1
        ii=i;
        for k=nz:-1:1
           c=mod(ii-1,nb_bin)+1; 
           colorIm(cptIm,k)=(c-0.501)*max_feat(k)./nb_bin; 
           ii=(ii-c)/nb_bin+1;
        end
%         c3=mod(i-1,nb_bin)+1;
%         ii=(i-c3)/nb_bin;
%         c2=mod(ii,nb_bin)+1;
%         c1=(ii-c2+1)/nb_bin+1;
        
%         colorIm(cptIm,1)=(c1-0.501)*256./nb_bin;
%         colorIm(cptIm,2)=(c2-0.501)*256./nb_bin;
%         colorIm(cptIm,3)=(c3-0.501)*256./nb_bin;
        cptIm=cptIm+1;
    end   
end

disp(['Creation of operators']);
%Creation of matrices H: pixels-> histograms
HlistePixel=1:1:nb_pixel;
HlisteOnes=double(ones(nb_pixel,1));

for i=1:nx
    for  j=1:ny
        index_pixel=nx*(j-1)+i;
        HlistePixel(index_pixel)=index_pixel;
        %Serch for which bin belong  I(c1,c2,c3)
        d=0;
        for k=1:nz
            d=d+(colorIm(:,k)-Ifeat(i,j,k)).^2;
        end
        
        [val index_color]=min(d);
        
        
        %Fast creattion of a sparse matrix with a list of index
        Hliste(index_pixel)=index_color;
        
        
    end
end

%Sparse operators H
H=sparse(Hliste,HlistePixel,HlisteOnes,nb_bin_totalIm,nb_pixel);

%Creation of operators A and B
op_H=@(u,k) H_ref{k}*sum(u(:));


%Creation of operators A^* and B^*
op_H_t=@(q,k) repmat(sum(H_ref{k}.*q(:)),[nb_pixel 1]);

%segmentation variable
u=double(ones(nx,ny,nb_label))/nb_label;

%Dual variables for histogram comparison
for k=1:nb_label
    p{k}=zeros(nb_bin_totalIm,1);
    q{k}=zeros(nb_bin_total_ref(k),1);
end

%Dual variables for TV
zx=zeros(nx,ny,nb_label);
zy=zeros(nx,ny,nb_label);

div=zeros(nx,ny,nb_label);


%% add a penalization inversely proportional to image gradient norm to fit image contours
%[gx gy]=gradient(I/255.);
%rho=repmat(rho./(1+0.1*sum(sqrt(gx.^2+gy.^2),3)),[1 1 nb_label]);


norme_H= (H*ones(nb_pixel,1));

for k=1:nb_label
    norme_Href{k}=nb_pixel*(H_ref{k}(:));%
    sigmaq{k}= (1.)./(norme_Href{k}+nb_bin_totalIm);
    sigmap{k}= (1.)./(norme_H+nb_bin_total_ref(k));
end

%Primal dual parameters
tau=1./(4.+2.);
sigmaz= 1./2.;

% initialize with local decision
figure(1);

u0=zeros(nx,ny);
for i=1:nb_pixel
    inde=find(H(:,i));
    maxi=H_ref{1}(inde);
    indi=1;
    for k=2:nb_label
        if H_ref{k}(inde)>maxi
            maxi=H_ref{k}(inde);
            indi=k;
        end
    end
    
    u0(i)=indi;
    
end


imagesc(double(u0)/nb_label);
colormap gray;
title('Initialization with highest local probability');
drawnow;


for k=1:nb_label
    Lline{k}=ones(2*nb_bin_totalIm*nb_bin_total_ref(k),1);
    Lcolumn{k}=zeros(2*nb_bin_totalIm*nb_bin_total_ref(k),1);
end

cpt_ref=ones(nb_label,1);
for i=1:nb_bin_totalIm
    for k=1:nb_label
        for j=1:nb_bin_total_ref(k)
            Lline{k}(cpt_ref(k))=(i-1)*nb_bin_total_ref(k)+j;
            Lline{k}(cpt_ref(k)+1)=Lline{k}(cpt_ref(k));
            cpt_ref(k)=cpt_ref(k)+2;
            
        end
    end
end


cpt_ref=ones(nb_label,1);
for k=1:nb_label
    for j=1:nb_bin_total_ref(k)
        for i=1:nb_bin_totalIm
            Lcolumn{k}(cpt_ref(k))=i;
            Lcolumn{k}(cpt_ref(k)+1)=nb_bin_totalIm+j;
            cpt_ref(k)=cpt_ref(k)+2;
        end
    end
end


for k=1:nb_label
    
    L{k}=sparse(Lline{k},Lcolumn{k},ones(2*nb_bin_totalIm*nb_bin_total_ref(k),1),nb_bin_totalIm*nb_bin_total_ref(k),nb_bin_totalIm+nb_bin_total_ref(k));
    r{k}=zeros(nb_bin_totalIm*nb_bin_total_ref(k),1);
    rt{k}=r{k};
end

ut=u;

disp(['Optimization']);
h = waitbar(0,'Initializing waitbar...');
for iter=1:niter_max
    
    waitbar(iter/niter_max,h,'In progress');
    if max(rho(:))>0
        zx(1:nx-1,:,:)= zx(1:nx-1,:,:)+sigmaz.*(ut(2:end,:,:)-ut(1:end-1,:,:));
        zy(:,1:ny-1,:)= zy(:,1:ny-1,:)+sigmaz.*(ut(:,2:end,:)-ut(:,1:end-1,:));
        
        
        normez=max(rho,sqrt(zx.^2+zy.^2));
        
        zx=rho.*zx./normez;
        zy=rho.*zy./normez;
    end
    
    for k=1:nb_label
        rr=L{k}'*rt{k};
        p{k}=p{k}+sigmap{k}.*(H*reshape(ut(:,:,k),[nx*ny,1]));
        q{k}=q{k}+sigmaq{k}.*(op_H(ut(:,:,k),k));
        p{k}=(sigmaq{k}.*p{k}-sigmap{k}.*q{k})./(sigmaq{k}+sigmap{k});
        p{k}(:)=p{k}(:)./(max(1.,abs(p{k}(:))));
        q{k}=-p{k};
        
    end
    
    
    
    if max(rho(:))>0
        div(1,:,:)=zx(1,:,:);
        div(2:end-1,:,:)=zx(2:end-1,:,:)-zx(1:end-2,:,:);
        div(end,:,:)=-zx(end-1,:,:);
        
        div(:,1,:)=div(:,1,:)+zy(:,1,:);
        div(:,2:end-1,:)=div(:,2:end-1,:)+zy(:,2:end-1,:)-zy(:,1:end-2,:);
        div(:,end,:)=div(:,end,:)-zy(:,end-1,:);
    end
    
    ut=u;
      
    
    for k=1:nb_label
        tmp=tau*(H'*p{k}+op_H_t(q{k},k)-reshape(div(:,:,k),[nx*ny 1]));
        u(:,:,k)=u(:,:,k)-reshape(tmp,[nx ny]);
    end
    %Projection onto simplex:
    if nb_label>2
        utmp=reshape(u, [nx*ny, nb_label]);
        utmp=simplex_projection(utmp);
        u=reshape(utmp,[nx ny nb_label]);
    else
        u(:,:,1)=max(0,min(1,(u(:,:,1)+1-u(:,:,2))/2));
        u(:,:,2)=1-u(:,:,1);
        
    end
    ut=2.*u-ut;
      
    if mod(iter,100)==0    
        [ a u0]=max(u,[],3);
        figure(1)
        imagesc(u0-1)
        
        colormap gray;
        title(['Segmentation map u']);
        drawnow;
    end
end

close(h)

% display results
colore{1}='r';
colore{2}='g';
colore{3}='b';
colore{4}='m';
figure;
u00=sum(u0,3);
imshow(I/255);
se=strel('square',2);
lw=3;

for k=nb_label:-1:1
    II=u00==k;
    hold on
    contour(imerode(II,se),colore{k},'Linewidth',lw);
end

