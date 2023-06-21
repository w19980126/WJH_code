%% ��ͬ����ͼ

ax1 = subplot(121);
set(ax1,'position',[0.01 0.05 0.42 0.95]);
imagesc(ax1,SPR_I);
axis off
axis square
set(gca,'fontsize',15,'fontweight','bold');
title({'��ɢ����λ'},'fontweight','bold');
caxis([0,0.1])
zlim([0 0.1])

ax2 = subplot(122);
set(ax2,'position',[0.45 0.05 0.42 0.95]);
imagesc(ax2,abs(Ir));
axis off
axis square
set(gca,'fontsize',15,'fontweight','bold');
title({'����ɢ�����Իָ�'},'fontweight','bold');

c = colorbar;
set(c,'position',[0.89 0.25 0.03 0.55])
caxis([0,0.1])

gFrame = getframe(gcf);
imwrite(gFrame.cdata,'C:\Users\20841\Desktop\tmp\�ĺ󼴷�.a.tif')
saveas(gcf,'C:\Users\20841\Desktop\tmp\�ĺ󼴷�.a.tif')
%% ����ͼƬ����
colormap(autumn)
set(gcf,'units','normalized');
Position = get(gcf,'position');
axis off
axis equal
set(gcf,'Position',Position);
set(gca,'fontsize',15,'fontweight','bold');
title({['k�ռ���ģʾ��ͼ']},'fontweight','bold');
xlabel('Horizontal coordinates(100nm)');
ylabel('Intensity');
box on
L = findobj(gca,'Type','Line');
set(L,'linewidth',1.5);
legend('�ο���λ','��ģȡ����','��ģ������')

colormap('parula');
xlim([0,41])
ylim([0,41])

%% ��������

savepath = 'F:\work\ScaterFeild\����\��תSPR\����\�����ǶȻָ�������';
savepath = fullfile(savepath,['k�ռ���ģʾ��ͼ'])
saveas(gcf,savepath)
figure
temp = phi;
temp(12:21) = -phi(12:21)+pi;
temp(22:31) = phi(22:31)-pi;
temp(32:end) = -phi(32:end)+2*pi;
plot(phi,temp);