%% Step 1: �����Բ�Ĳ��ڴ�ֱ���ϣ�����Ҫ��ת��ֵ��ʹ���8��

%% Step 2: ��������ͼ���ϻ��ƶ��ľ��룬������dist

%% Step 3: ��ֱ�����ƶ�ͼ����
%ʾ����������
siz=401;
im=peaks(siz);
IM=im+rot90(im,2);
dist=10*1;
temp=eye(siz);
moveMat=[temp(dist+1:end,:);temp(1:dist,:)];
IM2=moveMat*im+rot90(moveMat*im,2);
IM2=im*moveMat+rot90(im*moveMat,2);

%���
res = separateIM(IM,IM2,dist,'r');

%��ͼ
subplot(3,2,1)
imagesc(IM);title('��ת�ϲ�ͼ��1');
subplot(3,2,2)
imagesc(IM2);title('��ת�ϲ�ͼ��2');
subplot(3,2,3)
imagesc(im);title('ԭʼͼ��');
subplot(3,2,4)
imagesc(res);title('��ԭͼ��');
subplot(3,2,5)
imagesc((im-res)./im);title('ƫ��');

IM = IM1';
IM2 = IM2';
dist = 9;
%��ͼ
subplot(3,2,1)
imagesc(abs(IM));title('��ת�ϲ�ͼ��1');
subplot(3,2,2)
imagesc(abs(IM2));title('��ת�ϲ�ͼ��2');
subplot(3,2,3)
imagesc(abs(IM));title('ԭʼͼ��');
subplot(3,2,4)
imagesc(abs(res));title('��ԭͼ��');
subplot(3,2,5)
imagesc(abs(im-res)./abs(im));title('ƫ��');










