close all
clear all
clc
tiffpath = 'F:\work\散射场\实验数据\20220407_AuNPs_AgNWs\AuNPs_60nm\A1';
tiffs = dir(fullfile(tiffpath,'*.tiff'));

BG = MVMExtBG(tiffpath);
I = double(imread(fullfile(tiffpath,tiffs(2).name))) -BG;


[~,rectout] = imcrop(I/max(I,[],'all'));
rectout = round(rectout);
r1 = rectout(2);
r2 = rectout(2) + rectout(4);
c1 = rectout(1);
c2 = rectout(1) + rectout(3);
I = I(r1:r2,c1:c2);
close 

figure
imagesc(I);
axis off
axis equal
colormap(sunglow)

line([50 50+3000/73.8],[180 180],'color','white','linewidth',6)     % scalebar，长3um
line([65 225],[99 99],'linewidth',3,'color','white','linestyle','--')    % profile，中心左右各沿50个像素

F = fftshift(fft2(I));
figure
imagesc(abs(F))
axis off
axis equal
colormap(sunglow)

[center_raw,center_col,R,~] = findcircle(log(abs(F)),5,0,1);
peaks = [center_col,center_raw,R];

I1 = LineCut(I,[65 225],[99 99]);
close 
close
I1 = I1/max(I1);

temp = findobj(gca,'type','line');      % 在以往作图的基础上修改
temp.XData = 1:length(I1);
temp.YData = I1;

%% 理想低通
mask = EwaldMask(F,peaks,0.0001,'I',1);
I_flt = ifft2(ifftshift(F.*mask));
figure
imagesc(abs(I_flt))
axis off
axis equal
colormap(sunglow)

I1 = LineCut(abs(I_flt),[97 257],[123 123]);
close 
close 
I1 = I1/max(I1);

temp = findobj(gca,'type','line');      % 在以往作图的基础上修改
temp.YData = I1;

%% 理想低通
mask = EwaldMask(F,peaks,0.1,'I',1);
I_flt = ifft2(ifftshift(F.*mask));
figure
imagesc(abs(I_flt))
axis off
axis equal
colormap(sunglow)

I1 = LineCut(abs(I_flt),[65 225],[99 99]);
close 
close 
I1 = I1/max(I1);

temp = findobj(gca,'type','line');      % 在以往作图的基础上修改
temp.YData = I1;

%% 理想低通
mask = EwaldMask(F,peaks,0.1,'G',1);
I_flt = ifft2(ifftshift(F.*mask));
figure
imagesc(abs(I_flt))
axis off
axis equal
colormap(sunglow) 

I1 = LineCut(abs(I_flt),[97 257],[123 123]);
close 
close 
I1 = I1/max(I1);

temp = findobj(gca,'type','line');      % 在以往作图的基础上修改
temp.YData = I1;

%% 理想低通
mask = EwaldMask(F,peaks,0.2,'G',1);
I_flt = ifft2(ifftshift(F.*mask));
figure
imagesc(abs(I_flt))
axis off
axis equal
colormap(sunglow) 

I1 = LineCut(abs(I_flt),[97 257],[123 123]);
close 
close 
I1 = I1/max(I1);

temp = findobj(gca,'type','line');      % 在以往作图的基础上修改
temp.YData = I1;

%% 理想低通
mask = EwaldMask(F,peaks,0.26,'G',1);
I_flt = ifft2(ifftshift(F.*mask));
figure
imagesc(abs(I_flt))
axis off
axis equal
colormap(sunglow) 

I1 = LineCut(abs(I_flt),[97 257],[123 123]);
close 
close 
I1 = I1/max(I1);

temp = findobj(gca,'type','line');      % 在以往作图的基础上修改
temp.YData = I1;


%% 信噪比
for ii = 1:50
    f = ii/100;
    mask = EwaldMask(F,peaks,f,'G',1);
    I_flt = ifft2(ifftshift(F.*mask));
    S = LineCut(abs(I_flt),[97 257],[123 123]);
    
    S = S/max(S);
    N(1:76) = S(1:76);
    N(77:153) = S(84:end);
    SNR(ii) = S(77)/std(N);

    yy = S;   
    fun = @(par,xdata) par(1).*exp(-(xdata-par(2)).^2./(2*(par(3)).^2)) + par(4);   % Gaussian function
    par0 = [0.8,81,10,0];    % initial value
    x0 = 1:length(yy);
    par2 = lsqcurvefit(fun,par0,x0,yy');
    P(ii) = par2(3);

end
figure
plot(SNR)


yy = I1;   
fun = @(par,xdata) par(1).*exp(-(xdata-par(2)).^2./(2*(par(3)).^2)) + par(4);   % Gaussian function
par0 = [0.8,81,10,0];    % initial value
x0 = 1:length(yy);
par2 = lsqcurvefit(fun,par0,x0,yy');
figure
plot(x0,yy)
hold on
plot(x0,fun(par2,x0),'--');









