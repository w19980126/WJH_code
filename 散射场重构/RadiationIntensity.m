function [I,angleMask] = RadiationIntensity(M,peaks,scalefactor,varargin)
% 本函数用来得到k空间圆环的辐射图
% M：输入的k空间图片
% peaks：圆环中心，是一个一维向量[x,y]，其中x是圆心的横坐标，y是纵坐标
% scalefactor: 用来定义所截取圆环的粗细，在0~1之间
% varargin：选择性输入量，即角度分辨率，默认是30
% 角度是顺时针方向
    if isempty(varargin)
        N = 36;
    else 
        N = varargin{1};
    end

    siz = size(M);
    I = zeros(N,1);
    f0 = siz(2)/siz(1);

    mask = halfMask(siz,f0,peaks).*circlemaskF(siz,f0,peaks,scalefactor);
    M = abs(M).*mask;
    phi = angleFun(M,peaks).*mask;
    
    for ii = 1:N
        angleMask(:,:,ii) = (phi<2*pi/N*ii).*(phi>=2*pi/N*(ii-1)+eps);
        I(ii) = sum(squeeze(angleMask(:,:,ii)).*M,'all')/sum(angleMask(:,:,ii),'all');
    end
end

%% half mask
function mask = halfMask(siz,f0,peaks)
    center = (siz+1)/2;
    k0 = -1/((center(2)-peaks(1))/(f0*(center(1)-peaks(2)+eps)));   % 真实的k0
    [X,Y] = meshgrid(1:siz(2),1:siz(1));
    mask = sign(X - (k0*f0*(Y-center(1)) +center(2)));
    mask(mask>=0) = 0;
    mask(mask<0) = 1;
end

%% ideal lowpass filter
function mask = circlemaskF(siz,f0,peaks,scalefactor)
    D0 = scalefactor*peaks(3);
    [X,Y] = meshgrid(1:siz(2),1:siz(1));

    mask = sqrt((Y*f0-peaks(2)*f0).^2+(X-peaks(1)).^2) - peaks(3); 
    mask(abs(mask)<=D0) = 1;
    mask(abs(mask)>D0) = 0;
end

%% angle function
function phi = angleFun(M,peaks)
    siz = size(M);
    f0 = siz(2)/siz(1);

    [X,Y] = meshgrid(1:siz(2),1:siz(1));
    X = X+eps;Y = Y+eps;
    phi = atan(f0*(Y-peaks(2))./(X-peaks(1)));
    
    loc1 = logical((X<peaks(1)).* (Y>peaks(2)));
    phi(loc1) = pi + phi(loc1);
    loc2 = logical((X<peaks(1)).* (Y<peaks(2)));
    phi(loc2) = pi + phi(loc2);
    loc3 = logical((X>peaks(1)).* (Y<peaks(2)));
    phi(loc3) = 2*pi + phi(loc3);
end









