function fluid = initCompsADIFluid(props, varargin)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Initialize AD-solver fluid for compositonal model
%
% SYNOPSIS:
%   f = initCompsADIFluid()
%
% REQUIRED PARAMETERS:
%   relPerm - Relative Permeability
%   DENSITY-Density for water
%   PVTW    - Water PVT table
%   ROCK    - Rock table
%
% RETURNS:
% fluid - Fluid model include the function handle for calculating
% relative permeability, water density and rock compressibility
%}
        fluid = struct();
        % Phase equilibrium
        reg.SATNUM =[];
        reg.SATINX   =':';
        reg.PVTNUM =[];
        reg.PVTINX =':';
        reg.ROCKNUM = [];
        reg.ROCKINX = ':';
        % Permeability
        fns = fieldnames(props);
        for k = 1:numel(fns)
            fn = fns{k};
            if doAssign(fn)
                asgn = str2func(['assign',fn]);
                try
                    fluid = asgn(fluid, props.(fn), reg);
                catch ME
                    warning(msgid('Assign:Failed'), ...
                    'Could not assign property ''%s''. Encountered error: ''%s''',...
                        fn, ME.message);
                end
            end
        end

        function flag = doAssign(propNm)
        % Properties not resulting in individual functions
        excpt = {'SWL'   ,'SWCR'   ,'SWU' , ...
                 'SGL'   ,'SGCR'   ,'SGU' , ...
                 'SOWCR' ,'SOGCR'  , ...
                 'CNAMES','BIC'   , 'ACF', ...
                 'PCRIT' ,'TCRIT' , 'VCRIT',...
                 'MW',    'ZCRIT', ...
                 'ISWL'  ,'ISWCR'  ,'ISWU', ...
                 'ISGL'  ,'ISGCR'  ,'ISGU', ...
                 'ISOWCR','ISOGCR'};
        flag = ~any( strcmp(propNm , excpt) );
        end

end

