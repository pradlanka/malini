function [R]=hist2res2(H,fun)
% Evaluates Histogram data
% [R]=hist2res2(H)
%
% [y]=hist2res2(H,fun)
%	estimates fun-statistic
%
% fun	'mean'	mean
%	'std'	standard deviation
%	'var'	variance
%	'sem'	standard error of the mean
%	'rms'	root mean square
%	'meansq' mean of squares
%	'sum'	sum
%	'sumsq'	sum of squares
%	'CM#'	central moment of order #
%	'skewness' skewness 
%	'kurtosis' excess coefficient (Fisher kurtosis)
%
% see also: NaN/statistic
%
% REFERENCES:
% [1] C.L. Nikias and A.P. Petropulu "Higher-Order Spectra Analysis" Prentice Hall, 1993.
% [2] C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).
% [3] http://www.itl.nist.gov/
% [4] http://mathworld.wolfram.com/


%	Version 2.90
%	last revision 02 Apr 2002
%	Copyright (c) 1996-2002 by Alois Schloegl
%	e-mail: a.schloegl@ieee.org	

% .CHANGELOG
% 20.09.2001  calc of Quantiles improved, using FLIX.M    
% 16.02.2002  major revision, compatible to NaN/statistic 
% 01.03.2002  minor bug fix, works for HISTO3 and HISTO2 results 
% 02.03.2002  global FLAG_implicit_unbiased_estimation implemented
% 02.04.2002  bug fixed

if ~strcmp(H.datatype,'HISTOGRAM')
        fprintf(2,'ERROR: arg1 is not a histogram\n');
        return;
end;
if nargin<2, fun=[]; end;

global FLAG_implicit_unbiased_estimation; 
%%% check whether FLAG was already defined 
if exist('FLAG_implicit_unbiased_estimation')~=1,
	FLAG_implicit_unbiased_estimation=[];
end;
%%% set DEFAULT value of FLAG
if isempty(FLAG_implicit_unbiased_estimation),
	FLAG_implicit_unbiased_estimation=logical(1);
end;

sz 	= size(H.H)./size(H.X);
R.N 	= sumskipnan(H.H,1);
R.SUM 	= sumskipnan(H.H.*repmat(H.X,sz),1);
R.SSQ 	= sumskipnan(H.H.*repmat(H.X.*H.X,sz),1);
%R.S3P 	= sumskipnan(H.H.*repmat(H.X.^3,sz),1);	% sum of 3rd power
R.S4P 	= sumskipnan(H.H.*repmat(H.X.^4,sz),1);	% sum of 4th power
%R.S5P 	= sumskipnan(H.H.*repmat(H.X.^5,sz),1);	% sum of 5th power

R.MEAN	= R.SUM./R.N;
R.MSQ   = R.SSQ./R.N;
R.RMS	= sqrt(R.MSQ);
R.SSQ0  = R.SSQ-R.SUM.*R.MEAN;		% sum square of mean removed

if FLAG_implicit_unbiased_estimation,
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and STE are INF
else
    n1	= R.N;
end;

R.VAR  	= R.SSQ0./n1;	     		% variance (unbiased) 
R.STD  	= sqrt(R.VAR);		     	% standard deviation
R.SEM  	= sqrt(R.SSQ0./(R.N.*n1)); 	% standard error of the mean
R.SEV	= sqrt(n1.*(n1.*R.S4P./R.N+(R.N.^2-2*R.N+3).*(R.SSQ./R.N).^2)./(R.N.^3)); % standard error of the variance
R.Coefficient_of_variation = R.STD./R.MEAN;

R.CM2	= R.SSQ0./n1;
x       = repmat(H.X,sz) - repmat(R.MEAN,size(H.X,1),1);
R.CM3 	= sumskipnan(H.H.*(x.^3),1)./n1;
R.CM4 	= sumskipnan(H.H.*(x.^4),1)./n1;
%R.CM5 	= sumskipnan(H.H.*(x.^5),1)./n1;

R.SKEWNESS = R.CM3./(R.STD.^3);
R.KURTOSIS = R.CM4./(R.VAR.^2)-3;
R.MAD = sumskipnan(H.H.*abs(x),1)./R.N; % mean absolute deviation

H.PDF = H.H./H.N(ones(size(H.H,1),1),:);
R.ENTROPY = -sumskipnan(H.PDF.*log2(H.PDF),1);

if ~isempty(fun),
        fun = upper(fun);
        if strncmp(fun,'CM',2) 
                oo = str2double(fun(3:length(fun)));
                R = sumskipnan(H.PDF.*(x.^oo),1);
    	else	            
		R = getfield(R,fun);
	end;
end;

