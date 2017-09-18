function mu = LohrenzBrayClark(rhor, zi, T, Tc, Pc, MW)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Lohrenz-Bray-Clark correlation (1964) for gas and oil viscosities.
%
% SYNOPSIS:
%   mu = LohrenzBrayClark(rhor, zi, T, Tc, Pc, MW)
%
% REQUIRED PARAMETERS:
%   rhor - Reduced density, unitless
%     zi   - Fluid molar composition, ADI class
%     T   - Temperature, K
%     Tc  - Critical temperature, K
%     Pc  - Critical pressure, atm
%     MW-Molecular weight, g/mol
%
% RETURNS:
%     mu  -  Oil and gas viscosities, cp
%}
        Pc = Pc./atm;   
        a1 = 0.10230;
        a2 = 0.023364;
        a3 = 0.058533;
        a4 = -0.040758;
        a5 = 0.0093324;
     sumr1 = 0; sumr2 = 0; sumr3 = 0;
    for i = 1: length(zi)
            sumr1 =sumr1+ zi{i}.*Tc(i);  
            sumr2 = sumr2+zi{i}.*MW(i); 
            sumr3 =sumr3+ zi{i}.*Pc(i);
    end
    cosai = sumr1.^(1/6).*sumr2.^(-1/2).*sumr3.^(-2/3);

     sum1= 0; sum2 = 0;
    for i = 1: length(zi)
           Tri  = T./Tc(i); % Reduced temperature
            cosaii = Tc(i).^(1/6)./(MW(i).^0.5.*Pc(i).^(2/3));
            yetai = (double(Tri)<1.5).*34.*10^(-5).*(1./cosaii).*Tri.^0.94+ ...
            (double(Tri)>=1.5).*17.78.*10^-5.*(1./cosaii).*(4.58.*Tri-1.67).^(5./8);    
           sum1 = sum1+ (zi{i}.*yetai.*MW(i).^0.5);
           sum2 = sum2+ (zi{i}.*MW(i).^0.5);
    end
    yeta = sum1./sum2;
    mu = ((a1+a2.*rhor+a3.*rhor.^2+a4.*rhor.^3+a5.*rhor.^4 ).^4-10.^(-4))./cosai+yeta;
    mu  = mu.*milli; % Pa.s
end


