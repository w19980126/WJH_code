function R = TMM(modle,lambda,theta0,n,d)

    d = d/10^9;
    lambda = lambda/10^9;
    theta = zeros(length(n),1);
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
    R = (abs(r))^2;

end