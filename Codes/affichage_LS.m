addpath(genpath('../'))
load('Gold_Standards.mat');

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

n=2;
x=linspace(-5,5,100);
H=zeros(size(x)); H(x>0)=1;
y1=Heavyside_eta(x,1,2);
y2=Heavyside_eta(x,0.5,2);
y3=Heavyside_eta(x,0.1,2);
Y=[y3;y2;y1];
figure;
plot(x,H,'-b');
ylim([-0.5,1.5]);
xlabel("x"); ylabel("H(x)");
figure;
plot(x,Y)
legend('\epsilon = 0.1','\epsilon = 0.5','\epsilon = 1');
xlabel("x"); ylabel("H_{\epsilon}(x)");
