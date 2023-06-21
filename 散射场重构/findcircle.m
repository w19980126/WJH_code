function [center_raw,center_col,R,mask] = findcircle(I,thickness,fct,method)
% ������������ר�Ŷ�ȡk�ռ�Բ���ģ�����������������ʵ�ռ�ͼƬ
% ���ݱȲ�ͬ�������k�ռ���kx��kyƵ�ʿ̶ȵĲ���
% ������ͨ������Ƶ�ʹ�һ������f0����������ʵ�ռ�
% ͼƬ�ݺ�ȶ�����ĵ��ռ�Ƶ�ʿ̶Ȳ���
% ע�⣬��ʱ����İ뾶R����x������й�һ��֮��İ뾶������
% fcut��ָ�����ٰٷֱȵ����Ĳ��ֱ���ס����ֹ��Ƶ�ź�̫����ͼƬ�Աȶȼ�С

    [m,n] = size(I);
    f0 = n/m;   % Ƶ��У�����ӣ�ͨ��������ӽ�y����
                % Ƶ�ʷֱ��ʻ�ԭΪ��x����Ƶ�ʷֱ�����ͬ��ˮƽ
                % �������ж�Ϊ������ʱm<n����f0>1
    I(round(m/2-m*fct):round(m/2+m*fct),round(n/2-n*fct):round(n/2+n*fct)) = 0;
    figure
    imagesc(abs(I));
    axis off
    axis equal

    if nargin == 3 || (nargin == 4 && method == 0)     
        [x,y] = ginput();  % x�Ǹ����x���꣬y�Ǹ����y����
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
            pause(0.1);   % �۲��ֵ����Ľ���Ƿ�����
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
    y = y*f0;   % ��x����Ƶ�ʷֱ���Ϊ��׼����y�����Ƶ�ʷֱ��ʻ�ԭΪʵ��Ƶ�ʷֱ���
    [xc,yc,R] = MultiCircle(x,y);       % R����x����Ƶ�ʼ��Ϊ��׼
    center_raw = yc/f0;     % �˴���ԭΪ���ص�λ���Դ�Ϊ��ģ��Բ��
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
% �������ö�����Բ���
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
% ����������ͨ���õ���Բ�ĺͰ뾶�Լ��������ӻ��Բ����ģ���������ģһ��ûɶ��
    mask = zeros(size(I)); 
    for ii = 1:size(I,1)
        for jj = 1:size(I,2)
            r = sqrt((center_col - jj)^2 + (center_raw*f0 - ii*f0)^2);  % ����뾶ʱ�����ÿ͹���ʵ��Ƶ�ʣ�����Ҫ����f0
            if (r -R) <= thickness && (r - R) >= -thickness 
                mask(ii,jj) = 1;
            end
        end
    end
    
end