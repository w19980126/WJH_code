%% 不同的子图

ax1 = subplot(121);
set(ax1,'position',[0.01 0.05 0.42 0.95]);
imagesc(ax1,SPR_I);
axis off
axis square
set(gca,'fontsize',15,'fontweight','bold');
title({'点散射相位'},'fontweight','bold');
caxis([0,0.1])
zlim([0 0.1])

ax2 = subplot(122);
set(ax2,'position',[0.45 0.05 0.42 0.95]);
imagesc(ax2,abs(Ir));
axis off
axis square
set(gca,'fontsize',15,'fontweight','bold');
title({'忽略散射属性恢复'},'fontweight','bold');

c = colorbar;
set(c,'position',[0.89 0.25 0.03 0.55])
caxis([0,0.1])

gFrame = getframe(gcf);
imwrite(gFrame.cdata,'C:\Users\20841\Desktop\tmp\阅后即焚.a.tif')
saveas(gcf,'C:\Users\20841\Desktop\tmp\阅后即焚.a.tif')
%% 常规图片设置
colormap(autumn)
set(gcf,'units','normalized');
Position = get(gcf,'position');
axis off
axis equal
set(gcf,'Position',Position);
set(gca,'fontsize',15,'fontweight','bold');
title({['k空间掩模示意图']},'fontweight','bold');
xlabel('Horizontal coordinates(100nm)');
ylabel('Intensity');
box on
L = findobj(gca,'Type','Line');
set(L,'linewidth',1.5);
legend('参考相位','掩模取单环','掩模半置零')

colormap('parula');
xlim([0,41])
ylim([0,41])

%% 保存设置

savepath = 'F:\work\ScaterFeild\结论\旋转SPR\数据\两个角度恢复纳米线';
savepath = fullfile(savepath,['k空间掩模示意图'])
saveas(gcf,savepath)
figure
temp = phi;
temp(12:21) = -phi(12:21)+pi;
temp(22:31) = phi(22:31)-pi;
temp(32:end) = -phi(32:end)+2*pi;
plot(phi,temp);