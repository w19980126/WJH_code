% ---------------读入图片序列---------------------------
close all
clear all
tiffpath = 'D:\work\散射场\实验数据\20220413_Au_Colission\TIFF\B1';
tiffs = dir(fullfile(tiffpath,'*.tiff'));
expname = 'A1';

% ---------------获得mask---------------------------
N = 400;
tick = 1250;
% BG = MVMExtBG(tiffpath);
BG = double(imread(fullfile(tiffpath,tiffs(tick).name)));
temp = double(imread(fullfile(tiffpath,tiffs(N+tick).name))) - BG;

imagesc(temp)
% [x,y] = ginput(2);
% x = round(x);y = round(y);
% x = [min(x),max(x)];
% y = [min(y),max(y)];
% close gcf

x(1) = 1;x(2) = 639;
y(1) = 1;y(2) = 479;

F = fftshift(fft2(temp(y(1):y(2),x(1):x(2))));
[center_raw,center_col,R,~] = findcircle(log(abs(F)),5,0,1);
peaks = [center_col,center_raw,R];
mask = EwaldMask(F,peaks,0.1,'G',1);
clear F

% ---------------开始计算---------------------------
for ii = 1:N
    temp = double(imread(fullfile(tiffpath,tiffs(ii+tick).name))) - BG;
    I(:,:,ii) = temp(y(1):y(2),x(1):x(2));
    F(:,:,ii) = fftshift(fft2(temp(y(1):y(2),x(1):x(2))));
    Irec(:,:,ii) = abs(ifft2(ifftshift(mask.*squeeze(F(:,:,ii)))));
end

% ---------------演示---------------------------
figure
for ii = 1:N
    subplot(121)
    imagesc(squeeze(I(:,:,ii)))
    axis off
    axis equal
    title(num2str(ii))
    subplot(122)
    imagesc(squeeze(Irec(:,:,ii)));
    axis off
    axis equal
    title(num2str(ii))
    pause(0.01)
end
showSlide(Irec);
showSlide(I);

%% 保存为tiff文件
savepath = 'F:\work\散射场\实验数据\20220413_Au_Colission\Result\B1_firsttiff';
for ii = 1:N
    imwrite(uint16(squeeze(Irec(:,:,ii))),fullfile(savepath,[num2str(ii) '.tiff']));
end

%%
% ---------------保存为gif---------------------------
savepath = 'F:\work\散射场\实验数据\20220408_PsNPs_colission\PsNPs_300nm\Result\A1';
% mkdir(savepath)
gifroute = fullfile(savepath,['A1_1' '.gif']);
figure
for ii = 1:N
    subplot(121)
    imagesc(squeeze(I(:,:,ii)))
    axis off
    axis square
    caxis([-800 2500])
    subplot(122)
    imagesc(squeeze(Irec(:,:,ii)));
    axis off
    axis square
    caxis([200 1700])
    set(0,'defaultfigurecolor','w');
    Fgif = getframe(gcf);
    Igif = frame2im(Fgif);
    [Igif,map] = rgb2ind(Igif,256);
    if ii == 1
        imwrite(Igif,map,gifroute,'gif','Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(Igif,map,gifroute,'gif','DelayTime',0.1,'WriteMode','append');
    end
    pause(0.01)
end

% ---------------读取强度---------------------------
[int1,m] = gROI(I);
saveas(gcf,fullfile(savepath,[expname '_SPR']));
[int2] = gROI(Irec,m);
saveas(gcf,fullfile(savepath,[expname '_Reconstructioin']));


% ---------------保存为视频---------------------------
savepath = '‪H:\work\ScaterFeild\实验数据\20220322_colission\20220322_colission\TIFF\Result';
videoroute = fullfile(savepath,['A1' '.avi']);
myVideo = VideoWriter(videoroute);
open(myVideo)
figure
for ii = 1:N
    temp1 = squeeze(I(:,:,ii));
    temp1 = (temp1-min(temp1,[],'all'))/(max(temp1,[],'all')-min(temp1,[],'all'));
    temp1 = im2uint8(temp1);
    temp2 = squeeze(Irec(:,:,ii));
    temp2 = (temp2-min(temp2,[],'all'))/(max(temp2,[],'all')-min(temp2,[],'all'));
    temp2 = im2uint8(temp2);
    subplot(121)
    imshow(temp1)
    colormap(parula)
    caxis([0 255])
    subplot(122)
    imshow(temp2)
    colormap(parula)
    caxis([0 255])
    tempframe = getframe(gcf);
    tempframe = frame2im(tempframe);
    writeVideo(myVideo,tempframe);
end
close(myVideo);
  
% -----------------颗粒跟踪-------------------
figure
imagesc(squeeze(I(:,:,84)))
figure
imagesc(squeeze(I(:,:,84))>500)
temp = squeeze(I(:,:,100))>1800;
temp = medfilt2(temp);
figure
imagesc(temp)



temp1 = squeeze(I(:,:,84));
temp2 = imgaussfilt(temp1);
figure
imagesc(temp2)
temp3 = 1./(1+(abs(temp2)/500).^2.5);
figure
imagesc(temp3)
temp4 = temp3<0.7;
figure
imagesc(temp4)







