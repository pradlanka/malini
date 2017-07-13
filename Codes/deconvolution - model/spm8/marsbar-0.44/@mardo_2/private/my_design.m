function tf = my_design(des)
% returns 1 if design looks like it is of SPM99 type
% 
% $Id$
  
tf = 0;
if isfield(des, 'SPM'), des = des.SPM; end
if isfield(des, 'SPMid')
  tf = ~isempty(strmatch('SPM2', des.SPMid));
end