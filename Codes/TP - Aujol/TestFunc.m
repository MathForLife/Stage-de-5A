x=linspace(-5,5,100);
[X,Y]=meshgrid(x);

eta=0.1;
phi=1-sqrt(X.^2+Y.^2);

H1=Heavyside_eta(x,eta,1);
H2=Heavyside_eta(x,eta,2);
d1=delta_eta(x,eta,1);
d2=delta_eta(x,eta,2);

plot1D=false;
plot2D=false;

if plot1D==true
    figure();
    plot(x,H1,'-r',x,H2,'--b');
    figure();
    plot(x,d1,'-r',x,d2,'--b');
end

if plot2D==true
    H1=Heavyside_eta(phi,eta,1);
    H2=Heavyside_eta(phi,eta,2);
    d1=delta_eta(phi,eta,1);
    d2=delta_eta(phi,eta,2);

    subplot(221);
    contourf(X,Y,H1);
    title('H1');
    subplot(222);
    contourf(X,Y,H2);
    title('H2');
    subplot(223);
    contourf(X,Y,d1);
    title('d1');
    subplot(224);
    contourf(X,Y,d2);
    title('d2');
end

testMask=false;

if testMask==true
    Image=double(imread('eight.tif'));
    
    mask=roipoly(Image/255.);
    imshow(mask);

    phi=signed_distance_from_mask(mask);
    figure;
    surf(phi,'EdgeColor','none');
end

I=randn(10,10,10); z1=randn(10,10,10); z2=randn(10,10,10); z3=randn(10,10,10);
Ix=gradx(I); Iy=grady(I); Iz=gradz(I);
a=sum(sum(sum(Ix.*z1+Iy.*z2+Iz.*z3)))
b=sum(sum(sum(I.*-div(z1,z2,z3))))