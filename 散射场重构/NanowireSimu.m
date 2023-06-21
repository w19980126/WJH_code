close all
clear all
clc
%% 生成干涉散射场
theta = 70.5;
lambda = 670;
f_sca = SctCoef(lambda*1e-9);
% theta = rad2deg(f_sca.ResonanceAngle);
n = 1.5;
% kapa = f_sca.PropagationLength*1e6;         % 衰减长度，单位微米
kapa = 8;
phi = 0.1*pi;
scale_factor = 0.01;
M_size = 501;
theta_spp = 180;
theta_res = 70.5;
% theta_res = rad2deg(f_sca.ResonanceAngle);
[Ei,Es,F,I] = wave_generate(lambda,n,kapa,theta,phi,scale_factor,M_size,theta_spp,theta_res);

%% 生成纳米线
%--------------------------------------
% 纳米线
S = zeros(501);
sz = size(S);
alpha = 135;     % 纳米线沿ii方向倾斜角度
L = 45;    % 纳米线长度，单位像素
width = 0.5;      % 纳米线直径，单位像素
oriP = [251,251];   % 指定原点
for ii = 1:sz(1)
    for jj = 1:sz(2)
        x0 = (ii - oriP(1));
        y0 = (jj - oriP(2));
        k0 = tan(alpha/180*pi);
        d = abs(k0*x0 - y0)/(k0^2 + 1);
        if d <= width && sqrt(x0^2+y0^2)<L
            S(ii,jj) = 1;
        end
    end
end
% figure
% imshow(S)
% 
% N = 7;
% [x,y] = ginput(N);
% x = ceil(x);
% y = ceil(y);
% for ii = 1:N
%     S(y(ii),x(ii)) = 1;
% end
% figure
% imshow(S)

% -------------------------------------
% convolution
I1 = conv2(Es,S.*exp(i*angle(Ei)),'same');  % 复散射场
I2 = Ei;    % 平面参考波
I3 = conj(I1).*Ei + I1.*conj(Ei) + (abs(I1)).^2;    % 直流分量
I4 = conj(I1).*Ei + I1.*conj(Ei);   % 只考虑散射，忽略直流散射场项
% I4 = I4(101:400,101:400);
%% 重构并读取剖面强度

figure
imagesc(I4);
axis off
axis equal
colormap(sunglow)

F = fftshift(fft2(squeeze(I4)));
[center_raw,center_col,R,~] = findcircle(log(abs(F)),5,0,0);
peaks = [center_col,center_raw,R];

mask = EwaldMask(F,peaks,0.2,'G',1);
I_flt = ifft2(ifftshift(F.*mask));
figure
imagesc(2*abs(I_flt))       % 模值
axis off
axis equal
colormap(sunglow)

%% -----------------------------------
% 读取剖面强度
X = [200,300];Y = [200,300];
[I0,X,Y] = LineCut(I4,X,Y);
[I2,X,Y] = LineCut(2*abs(I_flt),X,Y);

hl = findobj(gca,'type','line');
hl(1).YData = I2/max(I2);
hl(2).YData = I0/max(I0);

LineCut(S,X,Y);

%% 曲线拟合
% -------------------------------------
% Gaussian function for profile fit 

I = I2;
yy = I2;   
yy = (yy - min(yy))/(max(yy) - min(yy));
fun = @(par,xdata) par(1).*exp(-(xdata-par(2)).^2./(2*(par(3)).^2)) + par(4);   % Gaussian function
par0 = [0.8,200,10,0];    % initial value
x0 = 1:length(yy);
par2 = lsqcurvefit(fun,par0,x0,yy');
figure
plot(x0,yy)
hold on
plot(x0,fun(par2,x0),'--');

% --------------------------------------
% Lorentzian function for profile fit
I = I2;
yy = I2;   % profile in y-direction
yy = (yy - min(yy))/(max(yy) - min(yy));
fun = @(par,xdata) par(1)./((par(2))^2+(xdata-par(3)).^2) + par(4);   % Gaussian function
par0 = [10,3,200,0];    % initial value
x0 = 1:length(yy);
par2 = lsqcurvefit(fun,par0,x0,yy');
figure
plot(x0,yy)
hold on
plot(x0,fun(par2,x0),'--');

%% 作图
figure
plot3(ones(size(I0)),1:length(I0),I2);
hold on
plot3(1.2*ones(size(I1)),1:length(I1),I1);
plot3(1.4*ones(size(I2)),1:length(I2),I0);















