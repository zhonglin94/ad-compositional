function [mix, names] = getCompsProps(zi, namesID)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Get properties for each component
%
% SYNOPSIS:
%   [mix, names] = getCompsProps(Z, namesID)
%
% REQUIRED PARAMETERS:
%      zi - Composition of the mixture
% namesID - Id for component in 'Components_db'
%
% RETURNS:
%     mix  - Mixture contains properties for each component
%  names - Name of the component
%}       
        propsRows = namesID;
        nc = length(propsRows);
        % Add path to database -- Do not change unless
        % File is moved from original position
        addpath(genpath('../compsDataBase')); 
        % Check if props database is loaded.
        if(~exist('Components_db','var')) 
        %Load as structured variable db.
        db = load ('Components_db.mat'); 
        %Transfer to cell array and clear db.
        props = db.Components_db; clear db; 
        end
        
        % Add binary interaction database
        if(~exist('BIC_db','var')) % Check if props database is loaded.
        db = load ('BIC_db.mat'); %Load as structured variable db.
        BIC = db.BIC_db; clear db; %Transfer to cell array and clear db.
        end
        
        % Binary interaction coefficient
        kijmat = zeros(nc,nc);
         for i = 1: nc
            for j = 1: nc
                kijmat(i,j) = BIC{propsRows(i),propsRows(j)};
            end
         end
        names = cell(1, nc);
        [Tc, Pc, w, MW, Vc] = deal(zeros(1,nc));
        for i = 1:nc
            %Error check to ensure proper withdrawal of compound i data
            err = 0;
            if (isempty(props{propsRows(i),2}))
                fprintf('Compound %d name not found.\n',i)
                err = 1;
            else
                %Read the name of compound i
                names{1,i} = props{propsRows(i),2};
            end
            if(isempty(props{propsRows(i),5}))
                fprintf('%s Critical Temperature not found.\n',namesA{i})
                err=1;
            else
                %Read the critical temperature of compound i
                Tc(i) = props{propsRows(i),5};
            end
            if(isempty(props{propsRows(i),4}))
                fprintf('%s Critical Pressure not found.\n', namesA{i})
                err=1;
            else
                %Read the critical pressure of compound i
                Pc(i) = props{propsRows(i),4};
            end
            if(isempty(props{propsRows(i),6}))
                fprintf('%s Acentric factor not found.\n', namesA{i})
                err=1;
            else
                %Read the Accentric factor of compound i
                w(i) = props{propsRows(i),6};
            end          
            if(isempty(props{propsRows(i),7}))
                fprintf('%s Mole weight not found.\n', namesA{i})
                err=1;
            else
                %Read the Critical volume of compound i
                MW(i) = props{propsRows(i),7};
            end        
            
            if(isempty(props{propsRows(i),10}))
                fprintf('%s Critical volume not found.\n', namesA{i})
                err=1;
            else
                %Read the crtical volume of compound i
                Vc(i) = props{propsRows(i),10};
            end         
            
            if(err == 1)
                disp('Terminating. Check props database row and folder path.')
                return;
            end
        end
        mix =   struct('Comp',zi,'Pc',Pc*atm,'Tc',Tc,'W',w,'MW',MW,'IntCoef',kijmat, 'Vc', Vc); 
end
