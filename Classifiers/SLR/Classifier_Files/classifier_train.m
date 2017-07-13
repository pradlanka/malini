function [Decision_Surf]=classifier_train(data,group_no,varargin)
if nargin > 2
    nlearn = varargin{1};
else
    nlearn = 300;
end
Decision_Surf = muclsfy_slrvarovrm_train(data, group_no,'nlearn',  nlearn, 'nstep', 100, ...
        'mean_mode', 'none', 'scale_mode', 'none', 'amax', 1e8);