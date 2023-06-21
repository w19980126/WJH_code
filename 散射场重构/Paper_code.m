
%% 生成干涉散射场
theta = 71;
n = 1.5;
lambda = 680;
kapa = 8;
phi = 0.1*pi;
scale_factor = 0.1;
M_size = 601;
theta_spp = 90;
theta_res = 70.5;
[Ei,Es,F,I] = wave_generate(lambda,n,kapa,theta,phi,scale_factor,M_size,theta_spp,theta_res);
imagesc(I)
figure
imagesc(abs(F))


%% Gaussian function for profile fit 

I = temp_I;
yx = abs(I(58,:));   % profile in x-direction
yy = abs(I(:,58));   % profile in y-direction
yx = (yx - min(yx))/(max(yx) - min(yx));    % normalization
yy = (yy - min(yy))/(max(yy) - min(yy));
fun = @(par,xdata) par(1).*exp(-(xdata-par(2)).^2./(2*(par(3)).^2)) + par(4);   % Gaussian function
par0 = [1000,58,10,0];    % initial value
x0 = 1:114;
par1 = lsqcurvefit(fun,par0,x0,yx);
par2 = lsqcurvefit(fun,par0,x0,yy');
figure
plot(x0,yx)
hold on
plot(x0,yy)
plot(x0,fun(par1,x0),'--');
plot(x0,fun(par2,x0),'--');

%% simulations for particles and nanowires

% to generate objects
S = zeros(501);

sz = size(S);
alpha = 135;     % 纳米线沿ii方向倾斜角度
L = 45;    % 纳米线长度，单位像素
width = 0.5;      % 纳米线直径，单位像素
oriP = [301,301];   % 指定原点
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
figure
imshow(S)

N = 5;
[x,y] = ginput(N);
x = ceil(x);
y = ceil(y);
for ii = 1:N
    S(y(ii),x(ii)) = 1;
end
figure
imshow(S)

% convolution
I1 = conv2(Es,S.*exp(i*angle(Ei)),'same');  % 复散射场
I2 = Ei;    % 平面参考波
I3 = conj(I1).*Ei + I1.*conj(Ei) + (abs(I1)).^2;    % 直流分量
I4 = conj(I1).*Ei + I1.*conj(Ei);   % 只考虑散射，忽略直流散射场项
figure
imagesc(I4)
figure
imagesc(I3)

F = fftshift(fft2(I4));
% [mask,peaks,deg]=GetFourierMask(log(abs(F)),10,0,-1);
[peaks(1,2),peaks(1,1),peaks(1,3),mask] = findcircle(log(abs(F)),10);
figure
imagesc(log(abs(F)))

mask = EwaldMask(F,peaks,0.9,1.2);
I_flt = (ifft2(ifftshift(F.*mask)));
figure
imagesc(2*real(I_flt))
figure
imagesc(2*abs(I_flt))

plot(S(272,:),'linewidth',1.5)
hold on
plot(2*real(I_flt(272,:)),'linewidth',1.5)
plot(2*abs(I_flt(272,:)),'linewidth',1.5)

axis off 
axis square
colormap(sunglow)

axis tight


%% 三张抠背景
tiffpath = 'F:\work\散射场\实验数据\20220407_AuNPs_AgNWs\AgNWs\A1';
tiffs = dir(fullfile(tiffpath,'*.tiff'));
NoiTiffs = zeros(3,480,640);
for ii = 1:3
    NoiTiffs(ii,:,:) = double(imread(fullfile(tiffpath,tiffs(ii).name)));
end
[BG,ref] = eliminateBg(NoiTiffs,1e-3,200);
I = squeeze(BG(2));
figure
imagesc(I)
F = fftshift(fft2(I));
imagesc(log(abs(F)));
% [mask,peaks,deg] = GetFourierMask(log(abs(F)),20,0,-1);
[center_raw,center_col,R,mask] = findcircle(log(abs(F)),5,0,0);
peaks = [center_col,center_raw,R];

% ----------------------------------
% 取两侧环内外roi
lineWidth = 10;
maskTemp = zeros(size(F));
temp=insertShape(maskTemp,'circle',peaks,'LineWidth',lineWidth,'Color','white');
mask = temp(:,:,1);
figure
imagesc(abs((ifft2((~mask.*F)))))

% -----------------------------------
% 对称取环内外roi
mask = EwaldMask(F,peaks,0.9,5);
mask(:,240:end) = 0;
I_flt = (ifft2(ifftshift(F.*mask)));
figure
imagesc(abs(I_flt))
axis off
axis square
figure
imagesc(mask)

% -----------------------------------
% 只取一侧环内roi
r = 0.75*peaks(1,3);
mask = zeros(size(squeeze(F)));
m2 = mask;
x = peaks(1,1);
y = peaks(1,2);
for ii = 1:size(F,1)
    for jj = 1:size(F,2)
        temp_r = sqrt((ii-y)^2+(jj-x)^2);
        if temp_r <= r
            mask(ii,jj) = 1;
        end
    end
end
figure
imagesc(abs((ifft2(mask.*F))))
axis off
axis square

figure
imagesc(I)
axis off
axis square
figure
imagesc(log(abs(mask.*F)))

% -------------------------------
% 手动取roi
temp = log(abs(F));
temp = (temp-min(temp,[],'all'))/(max(temp,[],'all')-min(temp,[],'all'));
mask = roipoly(temp1);
m2 = roipoly(temp1);
m3 = roipoly(temp1);
m4 = roipoly(temp1);
m = mask+m2+m3+m4;

figure
imagesc(abs(ifft2(ifftshift((m2).*F))));
axis off
axis square

%% 颗粒强度和重构结果与入射角度之间的关系
% ------------------------------------------------
% 读入图片
tiffpath = 'G:\work\ScaterFeild\实验数据\20220302_Colission\20220302_O3\TIFF\B_NPs_1';
tiffs = dir(fullfile(tiffpath,'*.tiff'));
N = round(length(tiffs)/2);
for ii = 1:N
    BG(:,:,ii) = double(imread(fullfile(tiffpath,tiffs(2*ii-1).name))); % 作为背景而采集的图像
    Figs(:,:,ii) = double(imread(fullfile(tiffpath,tiffs(2*ii).name))); % 作为信号而采集的图像
    I_sub(:,:,ii) = Figs(:,:,ii) - BG(:,:,ii);  % 信号减去背景后的图像
%     subplot(121)
%     imagesc(I_sub(:,:,ii));
%     axis off
%     axis equal
%     subplot(122)
%     imagesc(abs(fftshift(fft2(I_sub(:,1:480,ii)))));
%     axis off
%     axis equal
%     pause(0.5)
end

paper.NPs_SPR_Ang.BG = BG;
paper.NPs_SPR_Ang.Figs = Figs;
clear BG Figs 

% -----------------------------------------------
% 计算圆环圆心位置，并生成mask
figure
imagesc(I_sub(:,:,18));
[x,y] = ginput(2);
x = round(x);y = round(y);
Lcut = round(mean([abs(x(2)-x(1)),abs(y(2)-y(1))]));
Lcut = Lcut + mod(Lcut,2);
x = [min(x),min(x)+Lcut];
y = [min(y),min(y)+Lcut];
clear I_cut F_Ic mask
for ii = 1:N
    I_cut(:,:,ii) = I_sub(y(1):y(2),x(1):x(2),ii);
    F_Ic(:,:,ii) = fftshift(fft2(I_cut(:,:,ii)));
    [mask(:,:,ii),peaks(:,:,ii),deg(ii)] = GetFourierMask(log(abs(F_Ic(:,:,ii))),5,0,-1);
end

% -----------------------------------------------
% 取前35个数据对圆心位置进行线性插值
x1 = squeeze(peaks(1,1,1:35));
y1 = squeeze(peaks(1,2,1:35));
x2 = squeeze(peaks(2,1,1:35));
y2 = squeeze(peaks(2,2,1:35));
r = squeeze(peaks(1,3,1:35));
px1 = polyfit(1:35,x1,1);
px2 = polyfit(1:35,x2,1);
py1 = polyfit(1:35,y1,1);
py2 = polyfit(1:35,y2,1);
pr = polyfit(1:35,r,0);
plot(x1)
hold on
paper.NPs_SPR_Ang.peaks(1,1,1:N) = polyval(px1,1:N);
paper.NPs_SPR_Ang.peaks(1,2,1:N) = polyval(py1,1:N);
paper.NPs_SPR_Ang.peaks(2,1,1:N) = polyval(px2,1:N);
paper.NPs_SPR_Ang.peaks(2,2,1:N) = polyval(py2,1:N);
paper.NPs_SPR_Ang.peaks(1,3,1:N) = polyval(pr,1:N);
paper.NPs_SPR_Ang.peaks(2,3,1:N) = polyval(pr,1:N);

% -----------------------------------------------
% 生成圆形掩膜
r = 0.7*paper.NPs_SPR_Ang.peaks(2,3,1);
for kk = 1:N
    mask = zeros(size(squeeze(F_Ic(:,:,1))));
    x = paper.NPs_SPR_Ang.peaks(1,1,kk);
    y = paper.NPs_SPR_Ang.peaks(1,2,kk);
    for ii = 1:size(F_Ic,1)
        for jj = 1:size(F_Ic,2)
            temp_r = sqrt((ii-y)^2+(jj-x)^2);
            if temp_r <= r
                mask(ii,jj) = 1;
            end
        end
    end
    paper.NPs_SPR_Ang.MaskInsideEwald(:,:,kk) = mask;
end
paper.NPs_SPR_Ang.I_cut = I_cut;

% -----------------------------------------------
% 根据圆心位置画出掩膜，并基于此重构（掩膜只覆盖圆环，其余部分都保留）
figure
for ii = 1:N
    F_Ic = fftshift(fft2(paper.NPs_SPR_Ang.I_cut(:,:,ii)));
    maskTemp = log(abs(F_Ic));
    peaks = paper.NPs_SPR_Ang.peaks(:,:,ii);
    lineWidth = 20;
    temp=insertShape(maskTemp,'circle',peaks,'LineWidth',lineWidth,'Color','white');
    mask(:,:,ii) = ((temp(:,:,1) - maskTemp) ~= 0);
%     imagesc(abs(ifft2(ifftshift(~mask.*paper.NPs_SPR_Ang.F_Ic(:,:,ii)))))
%     axis off
%     axis equal
%     pause(0.5)
end

% -----------------------------------------------
% 对于每一幅图片都进行重构
figure
for ii = 1:N
    F_Ic = fftshift(fft2(paper.NPs_SPR_Ang.I_cut(:,:,ii)));
    paper.NPs_SPR_Ang.Recnst(:,:,ii) = abs(ifft2(ifftshift(F_Ic.*squeeze(paper.NPs_SPR_Ang.MaskInsideEwald(:,:,ii)))));
    imagesc(paper.NPs_SPR_Ang.Recnst(:,:,ii))
    title(num2str(ii))
    caxis([0 1200])
    colormap(sunglow)
    colorbar
    pause(0.1)
end

gROI(paper.NPs_SPR_Ang.Recnst);
% -----------------------------------------------
% 保存数据
filename = 'G:\work\ScaterFeild\实验数据\20220227AuNPs_Ag_NWs_angle\20220227\Result\A2.mat';
save(filename,'paper');

%% 撞击实验
% ------------------------------------------------
% 读入图片
close all
clear all
clc
tiffpath = 'F:\work\ScaterFeild\实验数据\20220227AuNPs_Ag_NWs_angle\20220227_2\TIFF\A7_nanowire';
tiffs = dir(fullfile(tiffpath,'*.tiff'));
N = round(length(tiffs))-1;
N = 200;
kk = 1;
for ii = 1:N
    I_sub(:,:,kk) = double(imread(fullfile(tiffpath,tiffs(kk+1).name))) ...
        - double(imread(fullfile(tiffpath,tiffs(1).name)));
    kk = kk + 1;
end

% -----------------------------------------------
% 计算圆环圆心位置，并生成mask
figure
imagesc(I_sub(:,:,N));
[x,y] = ginput(2);
x = round(x);y = round(y);
Lcut = round(mean([abs(x(2)-x(1)),abs(y(2)-y(1))]));
Lcut = Lcut + mod(Lcut,2);
x = [min(x),min(x)+Lcut];
y = [min(y),min(y)+Lcut];
close gcf
for ii = 1:N
    I_cut(:,:,ii) = I_sub(y(1):y(2),x(1):x(2),ii);
end
F = fftshift(fft2(squeeze(I_cut(:,:,N))));
% [mask,peaks,deg] = GetFourierMask(log(abs(F)),10,0,-1);
if size(peaks,1) == 2
    peaks = peaks(1,:);
end

[center_raw,center_col,R,mask] = findcircle(F,5);
peaks = [center_col,center_raw,R];

% ----------------------------------------------
% 取k空间一半的掩模
ii = 2;
F = fftshift(fft2(squeeze(I_cut(:,:,ii))));
mask = EwaldMask(F,peaks,0.8,1.2);
I_flt = (ifft2(ifftshift(F.*mask)));
figure
imagesc(abs(I_flt))
figure
imagesc(I_cut(:,:,ii))

figure
for ii = 1:N
    F_Ic = fftshift(fft2(paper.AuNP_Colission.I_cut(:,:,ii)));
    Recnst(:,:,ii) = abs(ifft2(ifftshift(F_Ic.*mask)));
    subplot(121)
    imagesc(Recnst(:,:,ii))
    title(num2str(ii))
    caxis([10 80])
    colormap(sunglow)
    colorbar
    subplot(122)
    imagesc(paper.AuNP_Colission.I_cut(:,:,ii))
    title(num2str(ii))
    caxis([-1000 1000])
    colormap(sunglow)
    colorbar
    pause(0.2)
end

% -----------------------------------------------
% 对于每一幅图片都进行重构
figure
for ii = 1:N
    F_Ic = fftshift(fft2(paper.AuNP_Colission.I_cut(:,:,ii)));
    paper.AuNP_Colission.Recnst(:,:,ii) = abs(ifft2(ifftshift(F_Ic.*squeeze(paper.AuNP_Colission.MaskInsideEwald))));
    subplot(121)
    imagesc(paper.AuNP_Colission.Recnst(:,:,ii))
    title(num2str(ii))
    caxis([10 700])
    colormap(sunglow)
    colorbar
    subplot(122)
    imagesc(paper.AuNP_Colission.I_cut(:,:,ii))
    title(num2str(ii))
    caxis([-2500 2500])
    colormap(sunglow)
    colorbar
    pause(0.1)
end

figure
for ii = 1:N
    subplot(121)
    imagesc(paper.AuNP_Colission.Recnst(:,:,ii))
    title(num2str(ii))
    caxis([10 500])
    colormap(sunglow)
    colorbar
    subplot(122)
    imagesc(paper.AuNP_Colission.I_cut(:,:,ii))
    title(num2str(ii))
    caxis([-2500 1500])
    colormap(sunglow)
    colorbar
    pause(0.1)
end

figure
for ii = 1:N
    F_Ic = fftshift(fft2(paper.AuNP_Colission.I_cut(:,:,ii)));
    paper.AuNP_Colission.Recnst(:,:,ii) = abs(ifft2(ifftshift(F_Ic.*squeeze(paper.AuNP_Colission.MaskInsideEwald))));
    y1 = mean(paper.AuNP_Colission.Recnst(240:280,:,ii));
    y2 = mean(paper.AuNP_Colission.I_cut(240:280,:,ii));
    yyaxis left
    plot(medfilt1(y1),'linewidth',2);
    ylim([50 1200]);
    yyaxis right
    plot(medfilt1(y2),'linewidth',2);
	ylim([-1800 800]);
    pause(0.2)
end

I1 = gROI(I_cut);
I2 = gROI(paper.AuNP_Colission.Recnst);
% -----------------------------------------------
% 保存数据
filename = 'G:\work\ScaterFeild\实验数据\20220302_Colission\20220302_O3\Result\D2.mat';
save(filename,'paper','-v7.3');

%% 纳米线数据处理
close all
clear all
clc
tiffpath = 'F:\work\ScaterFeild\实验数据\20220307_AgNWs\A7';
tiffs = dir(fullfile(tiffpath,'*.tiff'));

I1 = double(imread(fullfile(tiffpath,tiffs(1).name))) - double(imread(fullfile(tiffpath,tiffs(2).name)));
I2 = double(imread(fullfile(tiffpath,tiffs(2).name))) - double(imread(fullfile(tiffpath,tiffs(3).name)));

Icut1 = I1(1:479,1:479);
Icut2 = I2(1:479,1:479);
figure
imagesc(Icut1);
axis off
axis equal
colormap(sunglow)
figure
imagesc(Icut2);
axis off
axis equal
colormap(sunglow)
F = fftshift(fft2(squeeze(Icut1)));
[center_raw,center_col,R,mask] = findcircle(log(abs(F)),5);
peaks = [center_col,center_raw,R];


mask = EwaldMask(F,peaks,0.85,1.3);
I_flt = (ifft2(ifftshift(F.*mask)));
figure
imagesc(2*real(I_flt))
axis off
axis equal
colormap(sunglow)
figure
imagesc(2*abs(I_flt))
axis off
axis equal
colormap(sunglow)

[I0,X,Y] = LineCut(Icut1,X,Y);
[I1,X,Y] = LineCut(2*real(I_flt),X,Y);
[I2,X,Y] = LineCut(2*abs(I_flt),X,Y);
[I3,X,Y] = LineCut(Icut2,X,Y);

[I0,X,Y] = LineCut(Icut1);
[I2,X,Y] = LineCut(2*real(I_flt));
[I2,X,Y] = LineCut(2*abs(I_flt));
[I3,X,Y] = LineCut(Icut2);











