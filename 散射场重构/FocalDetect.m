%% 先采集一张图片，手动划取圆环，并得到最终的mask
tiffpath = 'F:\work\散射场\实验数据\20220315_60nmAuNPs_AgNWs\20220315_AgNWs_25fps\A1';
tiffs = dir(fullfile(tiffpath,'*.tiff'));

I = double(imread(fullfile(tiffpath,tiffs(5).name)));
[m,n] = size(I);
F = fftshift(fft2(squeeze(I)));
[center_raw,center_col,R,~] = findcircle(log(abs(F)),5);
peaks = [center_col,center_raw,R];

mask = EwaldMask(F,peaks,0.8,1.2);
I_flt = (ifft2(ifftshift(F.*mask)));
figure
imagesc(2*abs(I_flt))
axis off
axis equal

%% 利用上面得到的mask作为掩模实时处理传回的图片并显示


