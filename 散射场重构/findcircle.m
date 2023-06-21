function [center_raw,center_col,R,mask] = findcircle(I,thickness,fct,method)
% 本函数是用来专门读取k空间圆环的，本函数考虑了由于实空间图片
% 横纵比不同而引起的k空间中kx和ky频率刻度的差异
% 本函数通过引入频率归一化因子f0来消除由于实空间
% 图片纵横比而引入的倒空间频率刻度差异
% 注意，此时输出的半径R是以x方向进行归一化之后的半径！！！
% fcut是指定多少百分比的中心部分被盖住，防止低频信号太大导致图片对比度减小

    [m,n] = size(I);
    f0 = n/m;   % 频率校正因子，通过这个因子将y方向
                % 频率分辨率还原为与x方向频率分辨率相同的水平
                % 以行少列多为例，此时m<n，则f0>1
    I(round(m/2-m*fct):round(m/2+m*fct),round(n/2-n*fct):round(n/2+n*fct)) = 0;
    figure
    imagesc(abs(I));
    axis off
    axis equal

    if nargin == 3 || (nargin == 4 && method == 0)     
        [x,y] = ginput();  % x是各点的x坐标，y是各点的y坐标
    else
        while 1
            prompt = {'level:','block length:'};
            dlgtitle = 'Input';
            definput = {'0.98','5'};
            dims = [1 20];
            answer=inputdlg(prompt,dlgtitle,dims,definput);
            answer=str2double(answer);
            level=answer(1);
            blocklength=answer(2);
            Block=ones(size(I));
            Block(:,round(n/2-blocklength):end)=0;
            BW=ones(size(I));
            BW(I<level*max(I.*Block,[],'all')) = 0;
            BW=BW.*Block;
            
            F=figure;imagesc(BW);axis off;axis equal
            pause(0.1);   % 观察二值化后的结果是否满意
            answer=questdlg('Continue?','Continue?','YES','NO','YES');
            switch answer
                case 'YES'
                    close(F);
                    break;
                case 'NO'
                    close(F);
            end
        end
        [y,x]=find(BW==1);
    end
    hold on
    plot(x,y,'^')
    y = y*f0;   % 以x方向频率分辨率为基准，将y方向的频率分辨率还原为实际频率分辨率
    [xc,yc,R] = MultiCircle(x,y);       % R仍以x方向频率间隔为基准
    center_raw = yc/f0;     % 此处还原为像素单位，以此为掩模的圆心
    center_col = xc;
    mask = CircleMask(I,thickness,R,center_raw,center_col,f0);
    figure
    imagesc(abs(I).*mask)
    hold on
    plot(x,y/f0,'^r')
    axis off
    axis equal

end

%%
function [xc,yc,R] = MultiCircle(x,y)
% 本函数用多点进行圆拟合
    x_ = mean(x);
    y_ = mean(y);
    u = x-x_;
    v = y-y_;
    uuu = sum(u.^3);
    vvv = sum(v.^3);
    uu = sum(u.^2);
    vv = sum(v.^2);
    uv = sum(u.*v);
    uuv = sum(u.^2.*v);
    vvu = sum(v.^2.*u);
    uc = (uuv*uv -uuu*vv - vvu*vv + uv*vvv)/(2*(uv^2 - uu*vv));
    vc = (-uu*uuv + uuu*uv + uv*vvu - uu*vvv)/(2*(uv^2 - uu*vv));
    xc = uc + x_;
    yc = vc + y_;
    R = sqrt(mean((x-xc).^2 + (y-yc).^2));
end

%%
function mask = CircleMask(I,thickness,R,center_raw,center_col,f0)
% 本函数用来通过得到的圆心和半径以及缩放因子获得圆形掩模，但这个掩模一般没啥用
    mask = zeros(size(I)); 
    for ii = 1:size(I,1)
        for jj = 1:size(I,2)
            r = sqrt((center_col - jj)^2 + (center_raw*f0 - ii*f0)^2);  % 计算半径时考虑用客观真实的频率，故需要乘以f0
            if (r -R) <= thickness && (r - R) >= -thickness 
                mask(ii,jj) = 1;
            end
        end
    end
    
end