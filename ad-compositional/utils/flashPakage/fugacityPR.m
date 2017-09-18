function [fug, Z] = fugacityPR(comp, mix, p, t, phase)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Fugacity calculated by Peng-Robinson equation of state.
%
% SYNOPSIS:
%   [fug, Z] = fugacityPR(comp, mix, p, t, phase)
%
% REQUIRED PARAMETERS:
%   comp - Composition of the mixture, should be normalized to unity before inputed
%     mix - Properties of each inputted component
%     p   - Pressure of the mixture
%     t    - Temperature
%    phase- User defined. 0 for liquid while 1 for vapor
%
% RETURNS:
%     fug  - Fugacity of each component
%      Z   - Compressibility of liquid oil vapor
%}
        nc = numel(comp);
        if isnumeric(p)
            ng =numel(p);
        else
            ng = numel(p.val);
        end

        pci =   mix.Pc;      % Critical Pressures ( Pa )
        tci =   mix.Tc;          % Critical Temperatures ( K )
        wi  =   mix.W;         % Accentric Factors
        del =   mix.IntCoef; % Interaction Coefficients 

        R    = 8.314472;      %Pa.m^3/mol.K

        OmegaA  =   0.457235528921;
        OmegaB  =   0.0777960739039;

        aTc =   OmegaA.*R.^2.*tci.^2./pci;
        bTc =   OmegaB.*R.*tci./pci;
        % From 1979, Robinson, Robinson and Peng
        m   =   (wi <= 0.49).*(0.37464+1.54226*wi-0.26992*wi.^2) + ...
                (wi > 0.49).*(0.3796+1.485*wi-0.1644*wi.^2 + 0.01667*wi.^3);
        alpha=  ( 1 + m.* ( 1 - sqrt(t./tci) ) ).^2;

        b = 0; % Sum of b
        for i = 1 : nc
            b = b+ comp{i}.*bTc(:,i);
        end

        a = 0; % sum of a
        for i=1:nc
            for j=1:nc
                a=a+comp{i}.*comp{j}.*sqrt(aTc(i).*aTc(j).*alpha(i).*alpha(j) ).*(1-del(i,j));
            end
        end

        A=(a.*p)/(R.^2.*t.^2);
        B=(b.*p)/(R.*t);
        % Get the roots and find a appropriate one for the specified phase
        r = getCubicRoots(A,B);

        rootForSort = cell(1,3);
        for rootIndex = 1:3
               isImage = abs(imag(double(r{rootIndex})))>0;
               r{rootIndex}(isImage) = NaN;
               isNegative = double(r{rootIndex})<0;
               r{rootIndex}(isNegative) = NaN;         
               isNegative = double(r{rootIndex})<0;
               r{rootIndex}(isNegative) = NaN;  
               rootForSort{rootIndex} = double(r{rootIndex});        
        end

        [~,minIndex] = min(cell2mat(rootForSort),[],2);
        [~,maxIndex] = max(cell2mat(rootForSort),[],2);    
        % Phase = 0 --> GAS Phase,   phase = 1 --> LIQUID Phase
        Z = cell(ng,1);
        for i =1:ng
        Z{i}=phase.*r{minIndex(i)}(i)+(1-phase).*r{maxIndex(i)}(i);        
        end
        Z = vertcat(Z{:});

        psi = cell(1,nc);
        for ii=1:nc
            temp = 0;
            for jj=1:nc
                temp=temp+comp{jj}*sqrt(aTc(ii)*aTc(jj)*alpha(ii)*alpha(jj))*(1-del(ii,jj));
            end
            psi{ii} =temp;
        end
        fug = cell(1,nc);
        for ii=1:nc
            fug{ii} = comp{ii}.*p.*exp( bTc(ii)./b.*(Z-1)-log(Z-B)-A./(2.*sqrt(2).*B).*(2.*psi{ii}./a-bTc(ii)./b).*log((Z+2.414.*B)./(Z-0.414.*B)) );
        end        
end
   function roots = getCubicRoots(A,B)
        a =  1;
        b = -(1-B);
        c =  A-3.*B.^2-2.*B;
        d = -(A.*B-B.^2-B.^3);
        x1 = (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6*a.^2)).^2 + (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)... 
        - b./(3.*a) - (- b.^2./(9.*a.^2) + c./(3.*a))./(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6*a.^2)).^2 + (c./(3.*a) - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3)...
         - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3);
        x2 =  (- b.^2./(9.*a.^2) + c./(3.*a))./(2.*(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a) - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3)...
         - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)) + (3.^(1./2).*((- b.^2./(9.*a.^2) + c./(3.*a))./(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a)...
         - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3) + (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + ...
        (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)).*i)./2 - b./(3.*a) - ...
        (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)./2;

        x3 = (- b.^2./(9.*a.^2) + c./(3.*a))./(2.*(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a) - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3)... 
        - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)) - (3.^(1./2).*((- b.^2./(9.*a.^2) + c./(3.*a))./(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a)...
         - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3) + (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 +... 
        (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)).*i)./2 - b./(3.*a) - (((d./(2.*a) + b.^3./(27.*a.^3) ...
        - (b.*c)./(6.*a.^2)).^2 + (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)./2;
        roots = {x1,x2,x3};
   end
