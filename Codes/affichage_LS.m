addpath(genpath('../'))
load('Gold_Standards.mat');

coins_mask=Gold_Standards{2};
phi=signed_distance_from_mask(coins_mask);
surface(phi,'EdgeColor','none');
temp=abs(phi)==1;
hold on;
contour3(phi.*temp,'r');
hold off;
colorbar();
figure;
imshow(double(coins_mask));