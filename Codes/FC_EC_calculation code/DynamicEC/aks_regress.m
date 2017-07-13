% -----------------------------------------------------------------------
%   FUNCTION: aks_regress.m
%   PURPOSE:  perform multivariate regression
%
%   INPUT:  X           -   nvar (rows) by nobs (cols) observation matrix
%           NLAGS       -   number of lags to include in model
%
%   OUTPUT: ret.beta    -   coefficients 
%           ret.u       -   residuals 
%           ret.RSS     -   sum-square-error of residuals
%           ret.Z       -   covariance matrix of residuals
%        
%   Written AKS Sep 13 2004
%   Updated AKS December 2005
%   Ref: Seth, A.K. (2005) Network: Comp. Neural. Sys. 16(1):35-55
% -----------------------------------------------------------------------
function [ret] = aks_regress(X,nlags)

% figure regression parameters
[nvar,nobs] = size(X);

if(nvar>nobs) error('nvar>nobs, check input matrix'); end

% remove sample means if present (no constant terms in this regression)
m = mean(X');
if(abs(sum(m)) > 0.0001)
    mall = repmat(m',1,nobs);
    X = X-mall;
end


% construct lag matrices
lags = -999*ones(nvar,nobs-nlags,nlags);
for jj=1:nvar
    for ii=1:nlags
        lags(jj,:,nlags-ii+1) = X(jj,ii:nobs-nlags+ii-1);
    end
end

%  regression (no constant term)
regressors = [];
for ii=1:nvar
    regressors = [regressors squeeze(lags(ii,:,:))];
end
regressors = squeeze(regressors);
% bug fix Aug 8 2006 to allow nlags = 1
if size(regressors,1) == 1,
    regressors = regressors';
end

for ii=1:nvar
    xvec = X(ii,:)';
    xdep = xvec(nlags+1:end);
    beta(:,ii) = regressors\xdep;
    u(:,ii) = xdep-regressors*beta(:,ii);
    RSS(ii) = sum(u(:,ii).^2);
end

%   do r-squared
df_error = (nobs-nlags)-(nvar*nlags);
df_total = (nobs-nlags);
for ii = 1:nvar
    xvec = X(ii,nlags+1:end);
    rss2 = xvec*xvec';
    rss(ii) = 1 - (RSS(ii) ./ rss2);
    rss_adj(ii) = 1 - ((RSS(ii)/df_error) / (rss2/df_total) );
end

%   organize output structure
ret.beta = beta;
ret.u = u;
ret.rss = rss;
ret.rss_adj = rss_adj;
ret.Z = cov(u);

