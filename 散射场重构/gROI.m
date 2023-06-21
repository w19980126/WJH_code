function [I,varargout] = gROI(Tiffs,varargin)
% 从图片序列中划取ROI并返回强度曲线
% 必输入项是图片序列tiffs，必输出项是强度向量I
% 可选输入项是掩膜mask，可选输出项是掩膜mask
% 通过tempMat指定位置读取强度返回的mask没有意义
    
    [m,n,N] = size(Tiffs);
    
    if isempty(varargin)    % 若只输入图片序列，则需要人为地划定ROI，此时展示序列中的某一帧图片，在此图片中画ROI
        first_tiff = Tiffs(:,:,round(N/2));
        figure
        imagesc(first_tiff)
        axis off
        [x,y] = ginput(2);
        x = round(x);
        y = round(y);
        close gcf
        mask = zeros(size(first_tiff));
        mask(min(y):max(y),min(x):max(x)) = 1;
        
        I = linearRead(mask,Tiffs);
    else
        tempMat = varargin{1};      % varargin{1}是一个mask或者mask序列或者位置矩阵
        if size(tempMat,2) == n     % 输入一个mask或者mask序列
            mask = tempMat;
            I = linearRead(mask,Tiffs);
        else    % 否则表示输入的是一个位置序列pstn，记录每一帧中目标的实际位置，据此可以画出每一帧对应的mask
            r = varargin{2};
            cycle = varargin{3};
            [I,mask] = partMask(r,cycle,m,n,tempMat,Tiffs);
        end
    end
    
%     figure
%     plot(I,'linewidth',2);
    
    if nargout == 2
        varargout = {mask};
    end
     
end

%%
function I = linearRead(mask,Tiffs)
% 本函数是用来读取每张图片的强度
    
    if length(size(mask)) == 2
        I = zeros(size(Tiffs,3),1);
        N = size(Tiffs,3);
        for ii = 1:N
            temp_tiff = Tiffs(:,:,ii);
            I(ii) = sum(sum((mask.*temp_tiff)))/sum(sum(mask));
        end
    else
        I = sum(mask.*Tiffs,[1,2])./sum(mask,[1,2]);
        I = squeeze(I);
        I(isnan(I)) = 0;
    end
end

%%
function [I,mask] = partMask(r,cycle,m,n,tempMat,Tiffs)
% 本函数是用来读取通过tempMat指定的位置和帧数的颗粒的强度的，默认掩模半经是8个像素
    [X,Y] = meshgrid(1:n,1:m);      % X是列副本，Y是行副本
    frames = size(Tiffs,3);     % 每次处理多少帧
    mask = zeros(m,n,frames);
    ind1 = length(find(tempMat(:,4)<=(cycle-1)*frames));
    
    for ii = 1:frames
        ind1 = ind1 + ismember(ii+frames*(cycle-1),tempMat(:,4));   
        % 对应pstn矩阵中的第ind1行，若某帧未追踪到颗粒，则ind1跳过此帧不记录
        ind2 = ii + frames*(cycle-1);   
        % ind2对应的真实帧数，无论此时追踪到颗粒与否
        if ismember(ind2,tempMat(:,4))
            temp = sqrt((X-tempMat(ind1,2)).^2 + (Y - tempMat(ind1,3)).^2);
            temp(temp<=r) = 1;
            temp(temp>r) = 0;
        else
            temp = zeros(m,n);
        end
        mask(:,:,ii) = temp;
    end
    Tiffs(Tiffs<0) = 0;
    I = sum(mask.*Tiffs,[1 2])./sum(mask,[1 2]);    % 计算第一和第二维织成的面的均值
    I = squeeze(I);
    I(isnan(I)) = 0;
end























