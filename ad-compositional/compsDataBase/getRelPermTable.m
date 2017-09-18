function relPerm = getRelPermTable()
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%
% Get properties for each component
%
% SYNOPSIS:
%   relPerm = getRelPermTable()
%
% RETURNS:
%   relPerm - Relative permeability table
%} 

 addpath(genpath('../compsDataBase')); 
 % Check if props database is loaded.
  if(~exist('SWOF','var')) 
        %Load as structured variable db.
        db = load ('SWOF.mat'); 
        SWOF{1} = db.ans;
        clear db;
  end
  
   if(~exist('SGOF','var')) 
        db = load ('SGOF.mat'); 
        SGOF{1} = db.ans;
        clear db;
   end  
   
  assert(~isempty(SWOF), 'no SWOF is added');
  assert(~isempty(SGOF), 'no SGOF is added');  
  relPerm = struct('SWOF', SWOF,...
                                 'SGOF', SGOF);
end