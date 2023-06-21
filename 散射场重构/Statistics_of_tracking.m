%% 作图
%% 统计作图

spr = sprData;    % 未平滑的原始数据
rec = reconstructedData;  % 未平滑的原始数据
% for ii = 1:136
%     rec(:,ii) = rec(:,ii)/max(rec(:,ii));
%     spr(:,ii) = spr(:,ii)/max(spr(:,ii));
% end
n=50;
z=zeros(n+1,n+1);
z=sparse(z);
nonzrM = (spr ~= 0);	% 统计撞击数据中的非零元素none_zero_Matrix
dwell_time = sum(nonzrM,1);
spr(:,dwell_time>=50) = [];
rec(:,dwell_time>=50) = [];
% for ii = 1:size(rec,2)
%     rec(:,ii) = rec(:,ii)/max(rec(:,ii));
%     spr(:,ii) = spr(:,ii)/max(spr(:,ii));
% end
spr = nonzeros(spr);
rec = nonzeros(rec);
spr = spr/max(spr);
rec = rec/max(rec);
% maxspr = max(spr);
% maxrec = max(rec);
% minspr = min(spr(spr~=0),[],'all');
% minrec = min(rec(rec~=0),[],'all');
% rec = round(rec/maxrec*n);
% spr = round(spr/maxspr*n);
rec = round(rec*n);
spr = round(spr*n);

z1 = full(sparse(spr+1,rec+1,1,n+1,n+1));
z2 = full(sparse(rec+1,spr+1,1,n+1,n+1));

%%
% for ii = 1:size(spr,2)
% 
%     d1=spr(:,ii);
%     d2=rec(:,ii);
%     d1(d1==0)=[];
%     d2(d2==0)=[];
%     d1 = d1/spr(I1)*n;
%     d2 = d2/rec(I2)*n;
%     s = sparse(round(d1)+1,round(d2)+1,ones(size(d1)),n+1,n+1);
%     z = s+z;
% end
% 曲线拟合
[X,Y] = meshgrid(linspace(0,n+1,n+1),0:n);
w = z1;
A = sum(w.*X.^2,'all');
B = sum(w.*X,'all');
C = sum(w.*X.*Y,'all');
D = sum(w.*X,'all');
E = sum(w,'all');
F = sum(w.*Y,'all');
k = (C*E-B*F)/(A*E-B*D)
b = (A*F-C*D)/(A*E-B*D);

figure
hold on
% imagesc(log(z1*n/sum(z1,'all')))
imagesc(z1*n/sum(z1,'all'))
max(z1*n/sum(z1,'all'),[],'all')
colormap(violet)
fplot(@(x) k*x+b,[0 (n-b)/k],'r');
xlim([0 n])
ylim([0 n])

N = sum(z1,'all');
Y_bar = sum(w.*Y,'all')/N;
SYY = sum(w.*(Y-Y_bar).^2,'all')/N;
SSR = sum(w.*(Y-(k*X+b)).^2,'all')/N;
R2 = 1 - SSR/SYY;

%% 修改图片数据

h = get(gca,'Children');
him = h(2);
set(him,'CData',z/sum(z,'all'));
xlim([0 51])
axis square
colormap(violet)
caxis([0 1])
xlim([0 51])
ylim([0 51])