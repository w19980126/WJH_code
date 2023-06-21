classdef SctCoef
    %这个类用来计算散射场的参数而不输出参数
    %   properties
    %           wavelength          波长，单位m
    %           n                   各层折射率
    %           d                   各层厚度，单位m，要求与n长度和维度一致
    %           PropagationLength   传播长度，单位m
    %           ResonanceAngle      共振角，单位rad
    properties
        wavelength
        PropagationLength   % SPP传播长度，单位m
        ResonanceAngle      % 共振角，单位rad
        kspp
    end
    properties(Hidden,Constant)
        epsilon_au = -13.682 + 1.0356i;
        epsilon_water = 1.33^2;
        epsilon_glass = 1.55^2;
        c = 3e8;          
    end
    properties(Hidden)
        omega
        n 
        d
        k1
        k3      
    end
    
    methods
        function obj = SctCoef(wavelength)
            %UNTITLED 构造此类的实例
            %   此处显示详细说明
            obj.wavelength = wavelength;
            obj.omega = 2*pi*obj.c/wavelength;
            obj.k1 = 2*pi/wavelength*sqrt(obj.epsilon_glass);
            obj.k3 = 2*pi/wavelength*sqrt(obj.epsilon_water);
            obj.PropagationLength = ProLength(obj);
            obj.ResonanceAngle = ResAngle(obj);
            obj.kspp = Kspp(obj);
        end
        function ResonanceAngle = ResAngle(obj)
            obj.n = sqrt([obj.epsilon_glass obj.epsilon_au obj.epsilon_water]);
            obj.d = [1 50e-9 1];
            theta = linspace(0,pi/2,10000);
            R = zeros(size(theta));
            for ii = 1:length(R)
                R(ii) = TMMR('p',obj.wavelength*1e9,theta(ii),obj.d*1e9,obj.n);
            end
            [~,ind] = min(abs(R));
            ResonanceAngle = theta(ind);
        end
        function PropagationLength = ProLength(obj)
            im_kspp = imag(sqrt(obj.epsilon_water*obj.epsilon_au/(obj.epsilon_au + obj.epsilon_water) ...
                *obj.omega^2/obj.c^2));     % spp横向传播波矢虚部
            PropagationLength = 1/im_kspp;
        end
        function kspp = Kspp(obj)
            if isempty(obj.ResonanceAngle)
                obj.ResonanceAngle = ResAngle(obj);
                kspp = obj.k3*sin(obj.ResonanceAngle);
            else
                kspp = obj.k3*sin(obj.ResonanceAngle);
            end
        end
    end
end

