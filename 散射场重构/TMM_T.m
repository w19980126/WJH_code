function r = TMM_T(modle,lambda,theta0,n,d)

% 这是基于唐晋发的书编写的code，只能算反射率，透射率没看懂就没写
% modle，入射光偏振模式，modle有'TM'和'TE'两种
% lambda，入射光波长，单位是nm
% theta0，入射角度，单位是°
% n，各层折射率，形式是n = u + iv
% d，各层厚度，单位是nm

    d = d/10^9;
    lambda = lambda/10^9;
    theta = zeros(length(n),1);
    theta0 = theta0*pi/180;
    theta(1) = theta0;
    for ii = 2:length(n)
        theta(ii) = asin(n(1)/n(ii)*sin(theta(1)));
    end
    M = [1 0;0 1];
    temp_M = M;
    
    if strcmp(modle,'TM')
        for ii = 2:length(n)-1
            delta = 2*pi/lambda*n(ii)*d(ii)*cos(theta(ii));
            eta = n(ii)/cos(theta(ii));
            temp_M(1) = cos(delta);
            temp_M(2) = -1i*eta*sin(delta);
            temp_M(3) = -1i/eta*sin(delta);
            temp_M(4) = temp_M(1);
            M = M*temp_M;
        end
        eta_in = n(1)/cos(theta(1));
        eta_out = n(end)/cos(theta(end));
    elseif strcmp(modle,'TE')
        for ii = 2:length(n)-1
            delta = 2*pi/lambda*n(ii)*d(ii)*cos(theta(ii));
            eta = n(ii)*cos(theta(ii));
            temp_M(1) = cos(delta);
            temp_M(2) = -1i*eta*sin(delta);
            temp_M(3) = -1i/eta*sin(delta);
            temp_M(4) = temp_M(1);
            M = M*temp_M;
        end
        eta_in = n(1)*cos(theta(1));
        eta_out = n(end)*cos(theta(end));
    end
 
    M = M*[1;eta_out];
    Y = M(2)/M(1);
    r = (eta_in - Y)/(eta_in + Y);

end