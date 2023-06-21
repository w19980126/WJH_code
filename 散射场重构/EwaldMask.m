function mask = EwaldMask(F,peaks,scalefactor,modle,order,varargin)
% 此函数用来生成一个严格中心对称的mask，这个mask只遮蔽Ewald圆环
% F是k空间图像，其实只是用他的size而已
% scalefactor:缩放因子，用来衡量圆环展宽的，圆环展宽宽度是R*scalefactor，其中R是圆环半径

    if size(peaks,1) == 2
        peaks = peaks(1,:);
    end
    siz = size(F);
    f0 = siz(2)/siz(1);
    m_half = halfMask(siz,f0,peaks);

    if strcmp(modle,'I')
        mask_ring = ILPF(siz,f0,peaks,scalefactor,varargin);
    elseif strcmp(modle,'B')
        mask_ring = BLPF(siz,f0,peaks,scalefactor,order,varargin);
    elseif strcmp(modle,'G')
        mask_ring = GLPF(siz,f0,peaks,scalefactor,varargin);
    end
    mask = mask_ring.*m_half;

end

%% mask to cover half of the fig
function mask = halfMask(siz,f0,peaks)
    center = (siz+1)/2;
    k0 = -1/((center(2)-peaks(1))/(f0*(center(1)-peaks(2)+eps)));   % 真实的k0
    [X,Y] = meshgrid(1:siz(2),1:siz(1));
    mask = sign(X - (k0*f0*(Y-center(1)) +center(2)));
    mask(mask>=0) = 0;
    mask(mask<0) = 1;
end


%% ideal lowpass filter
function mask = ILPF(siz,f0,peaks,scalefactor,varargin)
    D0 = scalefactor*peaks(3);
    row = 1:siz(1);
    col = 1:siz(2);
    [X,Y] = meshgrid(col,row);
    D = hypot((X-peaks(1)),(Y-peaks(2))*f0) - peaks(3);    
    mask = ones(siz);

    mask(abs(D)<=D0) = 0;  
    if strcmp(varargin{1},'in')
        mask(D>0) = 0;
    end
end

%% Butterworth lowpass filter
function mask = BLPF(siz,f0,peaks,scalefactor,order,varargin)
    D0 = scalefactor*peaks(3);
    row = 1:siz(1);
    col = 1:siz(2);
    [X,Y] = meshgrid(col,row);
    D = hypot((X-peaks(1)),(Y-peaks(2))*f0) - peaks(3);    
    mask = 1 - 1./(1+(D/D0).^(2*order));
    if strcmp(varargin{1},'in')
        mask(D>0) = 0;
    end
end

%% Gauss lowpass filter
function mask = GLPF(siz,f0,peaks,scalefactot,varargin)
    D0 = scalefactot*peaks(3);
    row = 1:siz(1);
    col = 1:siz(2);
    [X,Y] = meshgrid(col,row);
    D = hypot((X-peaks(1)),(Y-peaks(2))*f0) - peaks(3);    
    mask = 1 - exp(-D.^2/(2*D0^2));
    if strcmp(varargin{1},'in')
        mask(D>0) = 0;
    end
end






