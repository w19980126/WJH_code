function [intensit,X,Y] = LineCut(I,pstx,psty)
% 本函数是用来得到图片I的强度剖面的
% 此函数支持手动取起点和终点，也支持输入起点和终点
% pstx指起点和终点的x坐标，是一个2*1的向量
% psty指起点和终点的y坐标，是一个2*1的向量
% 输出量intensity是强度的剖面
% X是起点和终点的x坐标
% Y是起点和终点的y坐标

    if nargin == 1
        figure
        imagesc(I)
        axis off
        axis equal
        colormap("gray")
        [x1,x2] = ginput(2);
        line(x1,x2,'Color','white','LineStyle','--','linewidth',3);
        colormap(gray)
        if nargout == 3
            X = x1;
            Y = x2;
        end
    elseif nargin == 3
        x1 = pstx;
        x2 = psty;
        if nargout == 3
            X = pstx;
            Y = psty;
        end
        figure
        imagesc(I)
        axis off
        axis equal
        colormap(gray)
        line(x1,x2,'Color','white','LineStyle','--','linewidth',3);
    end

    x1 = round(x1);     
    x2 = round(x2);
    
    if size(x1,1) == 1
        x1 = x1';
    end
    if size(x2,1) == 1
        x2 = x2';
    end

    temp = [x1,x2]';    
    x1 = temp(:,1);     % x1是第一个点的坐标，x1(1)是x列，x1(2)是y行
    x2 = temp(:,2);     % x2是第二个点的坐标
    [~,dim] = max(abs(x2 - x1));    % 确定切线的长度和方向性，以较长轴为参考进行延伸
    k = (x2(2)-x1(2))/(x2(1)-x1(1));
    n = max(abs(x2 - x1));  % 确定长轴方向与输出数据长度
    if k == inf % 划线与y轴平行
        raw1 = min([x1(2) x2(2)]);
        raw2 = max([x1(2) x2(2)]);
        col = x1(1);
        intensit = mean(I(raw1:raw2,(col-1):(col+1)),2);
    else
        intensit = zeros(n,1);     % 以较长轴进行延伸
        if dim == 1     % 划线更水平
            for ii = 1:n    
                loc = round([x1(1)+sign(x2(1)-x1(1))*ii;x1(2)+sign(x2(1)-x1(1))*k*ii]);         % x方向增加一个单位，y方向增加k个单位
                loc = sub2ind(size(I),loc(2),loc(1));
                intensit(ii) = mean(I([loc,loc+1,loc-1]));     % 三行进行平均
            end
        else    % 划线更倾斜
            for ii = 1:n
                loc = round([x1(1)+sign(x2(2)-x1(2))*ii/k;x1(2)+sign(x2(2)-x1(2))*ii]); % y方向增加一个单位，x方向增加1/k个单位
                intensit(ii) = mean([I(loc(2),loc(1)-1),I(loc(2),loc(1)),I(loc(2),loc(1)+1)]);   % 三列进行平均
            end
        end
    end
%     figure
%     plot(intensit,'linewidth',3);
%     xlabel('Pixel');
%     ylabel('Intensity');
%     title('剖面强度');
%     set(gca,'fontsize',15,'fontweight','bold');    
%     axis square

end