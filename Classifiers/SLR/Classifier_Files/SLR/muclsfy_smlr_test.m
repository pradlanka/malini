function [t_test_est] = muclsfy_slrvarovo_test(Decision_surf,x_test)
% Multi-class classification by sparse multinomial logistic regression (SMLR).
% 
% Normalization, parameter estimation, and performance evaluation are executed.
%
% -- Usage
% [ww, ix_eff_all, errTable_Tr, errTable_te, parm, AXall, Ptr, Pte] =...
%   muclsfy_smlr(x_train, t_train, x_test, t_test,  varargin)
%
% --- Input
% x_train :   [Nsamp_tr , Nfeat] 
% t_train :   [Nsamp_tr , 1]
% x_test  :   [Nsamp_te , Nfeat]
% t_test  :   [Nsamp_te , Nfeat]
%
% --- Optional Input
% parm = finputcheck(varargin, ...
%     {'scale_mode', 'string'  , {'all','each','stdall','stdeach','none'}, 'each';...
%      'mean_mode' , 'string'  , {'all','each','none'}, 'each';...
%      'ax0'       , 'real'    , []                   ,  [];...
%      'nlearn'    , 'integer' , [1 inf]              ,  1000;...
%      'nstep'     , 'integer' , [1 inf]              ,  100;...
%      'amax'      , 'real'    , [0 inf]              ,  1e8;...
%      'wmaxiter'  , 'integer' , [1 inf]              ,  50;...
%      'wdisp_mode', 'string'  , {'iter', 'off', 'final', 'notify'}, 'iter';...
%      'isplot'    , 'boolean' , []                   ,  0;...
%      'usebias'   , 'boolean' , []                   ,  1;...
%      'norm_sep'  , 'boolean' , []                   ,  0;... 
%      'displaytext','boolean' ,  []                  ,  1;...
%      'gamma0'    , 'real'    , [0 inf]              ,  0;...
%     });
%
% --- Output
% ww          :   Estimated weight parameters. [Nfeat, Nclass]
% ix_eff_all  :   Index of features survived. cell array of {Nclass}
% errTable_tr :   Counting table of each label estimated. [Nclass, Nclass]
% errTbale_te :   Counting table of each label estimated. [Nclass, Nclass]
% parm        :   Parmaters used in this routine. [struct]
% AXall       :   History of hyperparameters updating. [Nfeat*Nclass Nlearn]
% Ptr         :   Probaility of observing every label in training data. [Nsamp_tr Nclass]
%                 This value is used to put a label on each sample.
% Pte         :   Probaility of observing every label in training data. [Nsamp_te Nclass]
%                 This value is used to put a label on each sample.
%
% 2009/08/10 OY
% * Bug fix when 't_test' has only a single label of 't_train'.
% 2009/06/01 OY
% * The function is renamed (previously run_smlr.m)
%
% 2006/10/23 OY
% * 'Nclass' is removed from inputs.
% 2006/09/12 OY
%  * A field "nprobe" is introduced.
% 2006/09/06 OY  
%  * Output format modified (error table as output)
%  * Comment modified
% 2006/08/02 OY
% 2006/05/26 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

% char label -> number 
%[t_test] = label2num(t_test, label_names);

Nsamp_te = size(x_test,1);
%%% input check for optional parameter.
parm = Decision_surf.parm;
 
if ~isstruct(parm)
   error(parm);
end
parm.nsamp_te = Nsamp_te;
norm_sep = parm.norm_sep;
usebias  = parm.usebias;
scale = Decision_surf.scale;
base = Decision_surf.base;
ww =  Decision_surf.ww;

% normalize (sacling and baseline addjustment)
if norm_sep == 0
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode, scale, base);
else
    [x_test, scale, base] = normalize_feature(x_test, parm.scale_mode, parm.mean_mode);
end
% add a regressor for bias term
if usebias
    Xte = [x_test, ones(Nsamp_te,1)];
else
    Xte = [x_test];
end

%-----------------------
% Test
%----------------------
[t_test_est, ~] = calc_label(Xte, ww);
