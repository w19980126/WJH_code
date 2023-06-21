classdef Yu
    % 此类是用于复现余晖的工作的
    %   此处显示详细说明
    
    properties
        R                   % 以像素为单位的环半径
        I                   % 原始图片
        F                   % 原始图片的FFT
        kr                  % SPP波矢，是个矢量，故同时能确定其传播方向
        gamma               % 衰减系数，从实空间拟合得到，非从倒空间拟合
        Ep                  % 平面参考波
        h                   % 实空间点扩散函数
        M1                  % 掩膜一
        M2                  % 掩膜二
        Rcn                 % 重构结果
        FO                  % 滤波处理之后的k空间
    end
    
    methods
        function obj = Yu(I)
            obj.I = I;
        end
        function obj = k_gamma(obj,thickness,fct,method)
        % 获得波矢和衰减系数
            I_ = obj.I;
            F_ = fftshift(fft2(fftshift(I_)));
            [center_raw,center_col,R_,~] = findcircle(F_,thickness,fct,method);
            close;close;
            peaks = [center_col,center_raw,R_];
            xcenter = (size(I_,2)+1)/2;ycenter = (size(I_,1)+1)/2;
            f0 = size(F_,2)/size(F_,1);       % 校正因子，以x方向频率间隔为基准校正y方向的频率间隔
            [direc,~] = cart2pol(xcenter-peaks(1),f0*(ycenter-peaks(2)));
            kspp = 2*pi*R_/size(F_,2)*[cos(direc),sin(direc)];  % 以x方向频率间隔为准
            obj.kr = kspp*(kspp(1)<0) - kspp*(kspp(1)>=0);      % 默认是从右向左传播，否则校正之
            obj.gamma = Yu.ForGamma(I_);
            obj.F = F_;
            obj.R = R_;
        end
        function obj = ForMatrixes(obj,k1,ks1,ks2,varargin)
            % direct指定spp波的传播方向，direct默认为从右向左传播，strcmp(direct,'l')表示向左传播，
            % strcmp(direct,'r')表示向右传播，
            
            sz = size(obj.I);
            k_uint = 2*pi/sz(2);                                            % x方向频率间隔
            k1 = k1*k_uint;ks1 = ks1*k_uint;ks2 = ks2*k_uint;               % 输入的k值是像素单位的，这里转换一下
            
            if strcmp(varargin{1},'r')
                obj.kr = obj.kr*(obj.kr(1)<0) - obj.kr*(obj.kr(1)>=0);      % 默认是从右向左传播，否则校正之
            end
            X = -((sz(2)-1)/2):((sz(2)-1)/2);Y = -((sz(1)-1)/2):((sz(1)-1)/2);
            [X,Y] = meshgrid(X,Y);R_ = hypot(X,Y);
            obj.Ep = exp(1i*(obj.kr(1)*X + obj.kr(2)*Y));            % 参考波
            
            % ------------------ 计算点扩散函数及相关参数 --------------
            obj.h = exp((1i*vecnorm(obj.kr)-1/obj.gamma)*R_);        % 实空间点扩散函数
            Fh = fftshift(fft2(fftshift(obj.h)));                    % 频域传递函数，加eps防止除以极小数
            Fh(Fh == 0) = eps;                                       % 不知道为什么即便加上eps也总有0值，此处为止强行赋值
            
            % ------------------ 计算M1 -----------------------
            kx = 2*pi*X/sz(2);ky = 2*pi*Y/sz(1);k0 = vecnorm(obj.kr);
            obj.M1 = exp(-power((hypot(kx,ky)-k0)/ks1,2));     
            
            % ------------------ 计算M2 -----------------------
            if strcmp(varargin,'r')
                tempM = exp(-power(kx-k1,2)/ks2^2);
                tempM(kx<k1) = 1;
            else
                tempM = exp(-power(kx+k1,2)/ks2^2);
                tempM(kx>-k1) = 1;
            end
            obj.M2 = tempM;
            F_Es = fftshift(fft2(fftshift(obj.I.*obj.Ep)));
            FO_ = F_Es./Fh.*obj.M1.*obj.M2;
            obj.FO = FO_ + eps*(1+1i);
            obj.Rcn = abs(ifftshift(ifft2(ifftshift(FO_))));
        end
    end
    
    methods(Static)
        function gamma = ForGamma(I)
            % 从中心向尾巴取点
            [intensity,X,Y] = LineCut(I);
            close gcf;
            [~,R] = cart2pol(X(2)-X(1),Y(2)-Y(1));
            StartInd = find(intensity == max(intensity));
            x = linspace(0,R*(length(intensity) - StartInd)/(length(intensity)-1),(length(intensity) - StartInd+1))';
            y = intensity(StartInd:end);y = y/max(y);
            p0 = [1,length(x),0,0];
            f = @(p,x) p(1)*exp(-(x-p(4))/p(2)) + p(3);    
            p = lsqcurvefit(f,p0,x,y);
            figure
            plot(x,y)
            hold on
            plot(x,f(p,x)) 
            gamma = p(2);        
        end
    end
end

