%% 生成干涉散射场
clear;
theta = 70.5;
n = 1.5;
lambda = 680;
kapa = 8;
phi = pi;
scale_factor = 0.1;
M_size = 501;
theta_spp = 0;
theta_res = 70.5;
[Ei,Es,F,I] = wave_generate(lambda,n,kapa,theta,phi,scale_factor,M_size,theta_spp,theta_res);
figure;imagesc(I);axis off; axis equal; colormap('gray'); 

%% 重构
[center_raw,center_col,R,~] = findcircle((abs(F)),5,0,0);
peaks = [center_col,center_raw,R];

% mask = EwaldMask(F,peaks,0.3,'G',1,'in');
mask = EwaldMask(F,peaks,0.2,'G',1);

I_rcn = ifftshift(ifft2(ifftshift(F.*mask)));
figure
imagesc(abs(I_rcn))
axis off
axis equal
colormap(sunglow)

ax2 = axes('position',[0.555 0.575 0.35 0.35]);
imagesc(abs(I_rcn(220:282,220:282)))
axis square
axis off    % 调完位置再关掉坐标轴
