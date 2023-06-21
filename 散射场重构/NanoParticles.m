%% 纳米颗粒数据处理
%% 中值滤波法去除背景，无明场对比
close all
clear all
clc
tiffpath = 'D:\work\散射场\实验数据\20220408_PsNPs_colission\PsNPs_300nm\BG';
tiffs = dir(fullfile(tiffpath,'*.tiff'));

BG = MVMExtBG(tiffpath);

I = double(imread(fullfile(tiffpath,tiffs(1).name))) - BG;

figure
imagesc(I);
axis off
axis equal
colormap(sunglow)

F = fftshift(fft2(I));
[center_raw,center_col,R,~] = findcircle(log(abs(F)),5,0,0);
peaks = [center_col,center_raw,R];

mask = EwaldMask(F,peaks,0.15,'G',2);
I_rcn = abs(ifft2(ifftshift(F.*mask)));
figure
imagesc(I_rcn)
axis off
axis equal
colormap(sunglow)

line([50 50+3000/73.8],[250 250],'color','white','linewidth',6)     % scalebar，长5um


%% 读取剖面强度

[int0,X,Y] = LineCut(I,X,Y);
[int1,X,Y] = LineCut(2*abs(I_rcn),X,Y);
[int1,X,Y] = LineCut(2*real(I_rcn),X,Y);
[int1,X,Y] = LineCut(2*imag(I_rcn),X,Y);
[int2,X,Y] = LineCut(I2,X,Y);

[int0,X,Y] = LineCut(I);
[int1,X,Y] = LineCut(2*abs(I_rcn));
[int2,X,Y] = LineCut(I2);

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
yy = I2(100:150);   % profile in y-direction
yy = (yy - min(yy))/(max(yy) - min(yy));
figure
plot(yy)
fun = @(par,xdata) par(1)./((par(2))^2+(xdata-par(3)).^2) + par(4);   % Gaussian function
par0 = [10,7,32,0];    % initial value
x0 = 1:length(yy);
par2 = lsqcurvefit(fun,par0,x0,yy');
figure
plot(x0,yy)
hold on
plot(x0,fun(par2,x0),'--');

syms x
eqn = fun(par2,x) == 0.1*fun(par2,par2(3));
s = double(solve(eqn,x));
78*(s(2)-s(1))

%% 曲线作图
temp0 = I0/(max(I0));
temp1 = I1/(abs(min(I1)));
temp2 = I2/(max(I2));
figure
plot3(ones(size(I0)),1:length(I0),temp2);
hold on
plot3(1.22*ones(size(I1)),1:length(I1),temp1);
plot3(1.4*ones(size(I2)),1:length(I2),temp0);

%% 中值滤波法去除背景，无明场对比
close all
clear all
clc
tiffpath = 'F:\work\散射场\实验数据\20220408_PsNPs_colission\PsNPs_300nm\BG';
tiffs = dir(fullfile(tiffpath,'*.tiff'));

BG = MVMExtBG(tiffpath);

I = double(imread(fullfile(tiffpath,tiffs(2).name))) - BG;

[~,rectout] = imcrop(I/max(I,[],'all'));
rectout = round(rectout);
r1 = rectout(2);
r2 = rectout(2) + rectout(4) +(mod(rectout(4),2) ~= 0);
c1 = rectout(1);
c2 = rectout(1) + rectout(3)+(mod(rectout(3),2) ~= 0);
I = I(r1:r2,c1:c2);
close 

figure
imagesc(I);
axis off
axis equal
colormap(sunglow)

F = fftshift(fft2(I));
[center_raw,center_col,R,~] = findcircle(log(abs(F)),5,0,1);
peaks = [center_col,center_raw,R];

mask = EwaldMask(F,peaks,0.15,'G',2);
I_rcn = (ifft2(ifftshift(F.*mask)));
figure
imagesc(abs(I_rcn))
axis off
axis equal
colormap(sunglow)

line([40 40+2000/73.8],[165 165],'color','white','linewidth',6)     % scalebar，长5um

%% 利用上一节得到的mask和BG，自动处理图片并保存为600dip的tiff格式
savepath = 'F:\work\散射场\实验数据\20220407_AuNPs_AgNWs\AgNWs\Result\A1_减去第一帧';
mkdir(savepath)
hwait = waitbar(0);
for ii = 1:length(tiffs)
    I = double(imread(fullfile(tiffpath,tiffs(ii).name))) - BG;
    F = fftshift(fft2(I(y(1):y(2),x(1):x(2))));

    f1 = figure('visible','off');
    imagesc(I(y(1):y(2),x(1):x(2)))
    axis off
    axis normal
    colormap(sunglow)
    print(f1,fullfile(savepath,[num2str(ii) '_SPR']),'-dtiff','-r600');
    close(f1)

    f2 = figure('visible','off');
    I_rcn = ifft2(ifftshift(F.*mask));
    imagesc(2*abs(I_rcn))
    axis off
    axis normal
    colormap(sunglow)
    print(f2,fullfile(savepath,[num2str(ii) '_重构']),'-dtiff','-r600');
    close(f2)

    waitbar(ii/length(tiffs),hwait,[num2str(ii) '/' num2str(length(tiffs))]);
end
delete(hwait)




