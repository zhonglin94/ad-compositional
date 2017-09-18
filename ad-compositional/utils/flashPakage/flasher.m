function [ ki, alpha, st, Z] = flasher(model, zi , p, varargin)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Substitution method for flash calculation and negative flash for phase identification
%
% SYNOPSIS:
%   [ ki, x0, st] = flasher(model, zi , p, varargin)
%
% REQUIRED PARAMETERS:
% model - Compositional model
%     zi - Composition of the mixture
%     p - Pressure of the mixture
%
% OPTIONAL PARAMETERS
%     t  - User specified temperature, Kelvin. By defaut, t = tres
%
% RETURNS:
%     ki  -  Equilibrium constant at equilibrium
%    alpha-Vapor to liquid ratio at equilibrium
%     st  - Status for the mixture
%            only oil, ki = NaN; alpha = 0;  st =1
%            only gas, ki =NaN; alpha = 1;  st =2
%            gas and oil, ki = ki; alpha = alpha; st =3
%}
          t = model.inputdata.t ; 
         if numel(varargin) >= 1
             t = varargin{1}; 
         elseif numel(varargin) >= 2
             disp('flasher: Only temperature can be specified, you may inputed too many arguments');
         end
        mix = model.inputdata.props.mix;
        pci =   mix.Pc;         % Critical Pressures ( Pa )
        tci =   mix.Tc;          % Critical Temperatures ( K )
        wi  =   mix.W;         % Accentric Factors
         % Use Wilson equation to initialize vapor/liquid equilbrium constant
        ki  =   (pci./p).*exp( 5.37*(1+wi).*(1-tci./t) );
        x0  =   0.5;      % Initial guess for V/F ratio

        es  =   1e-12;    % Error for Newton Method
        n   =   100;       % Maximum number of iterations
        flag=   true;     % Entering the while loop   

        fV=ones(size(zi));
        fL=ones(size(zi));

        while (1/length(zi)*sum(log(fV./fL).^2) > 1e-12 || flag)
            flag=false;
            VFRatio
            alpha=x0;
            % Calculating Vapor and Liquid phase compositions
            xi  =   zi./(1+(ki-1)*alpha);
            yi  =   ki.*xi;
            % Normalizing Values
            xi  =   xi/sum(xi);
            yi  =   yi/sum(yi);
            xi = num2cell(xi,1);
            yi = num2cell(yi,1);
            % 0 for gas; 1 for oil
        %     [fL  ZL]=model.fugacity (xi,p,t,1);
        %     [fV ZV]=model.fugacity(yi,p,t, 0);
             [fL  ZL]=model.fugacity (xi,p,t,1);
             [fV ZV]=model.fugacity(yi,p,t, 0);
             fL = cell2mat(fL);fV = cell2mat(fV);
            % Update K-values
            ki=ki.*fL./fV;
        end
        % Reporting values
        if ~isreal(alpha)
         error('flasher: calculate vapor liquid ratio is unreal, maybe due to unfavorable pressure or compositions');
        elseif alpha<0
                st = 1; % Liquid like
                zi = num2cell(zi,1);
                [~, ZL]=model.fugacity(zi, p,t,1);
                % We them the same compressibility factor even though gas is not present
                Z = [ZL ZL]; 
                alpha = eps;
                ki = ones(1, numel(zi)).*eps; % Avoid NaN value
         elseif alpha>1
                st = 2; % Vapor like
                zi = num2cell(zi,1);
                [~, ZV]=model.fugacity(zi, p,t,1);
                % We them the same compressibility factor even though oil is not present
                Z = [ZV, ZV];
                alpha = 1;
                 ki = ones(1, numel(zi))./eps;
        else 
               st = 3;  % Vapor liquid like
               Z = [ZL, ZV];
        end

    function VFRatio
        % Determing V/F ratio
        iter=0;
        persistent io
        persistent q
        if isempty(io)
            io=0;
        end
        if isempty(q)
            q=0;
        end
        fx=sum( ((ki-1).*zi)./(1+(ki-1)*x0) );
        dfx=-1*sum( ((ki-1).^2.*zi)./((1+(ki-1)*x0).^2));
        ea=es+eps;
        while (ea>es)&&(iter<=n)&&(dfx~=0)
            x1=x0-fx/dfx;
            err=abs(x1-x0);
            ea=abs((x1-x0)/x1)*100;
            x0=x1;
            fx=sum( ((ki-1).*zi)./(1+(ki-1)*x0) );
            dfx=-1*sum( ((ki-1).^2.*zi)./((1+(ki-1)*x0).^2));
            iter=iter+1;
        end
        q=q+1;
    end
end
