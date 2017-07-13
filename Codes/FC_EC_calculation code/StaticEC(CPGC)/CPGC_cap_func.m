
function [con ins] = CPGC_cap_func(Y,X,lags)


if nargin==3, k=0; Ext=[];
elseif nargin==4, k=leads; Ext=[];
end
q=lags;
if size(Y,1)~=size(X,1), 
	error('The regressands and regressors must have the same sample size');
end
[T p]=size(X);
[T n]=size(Y);

XX=[];
if k==-1
	Y(1:q,:)=[];
	for j=1:q,
		XX=[X(q+1-j:T-j,:) XX];
	end
else
	Y([1:q T+1-k:T],:)=[];
	for j=1:k+q+1
		XX=[X(q+k+2-j:T+1-j,:) XX]; % Columns ordered from lags to leads
    end
end
[TT M]=size(XX);
%XX=[ones(TT,1) XX]; % Constant term included in the 1st column
%[TT M]=size(XX);

% Include additional variables, if specified, into the last columns of XX
if ~isempty(Ext),
	if size(Ext,1)~=size(XX,1),
		error('Ext must have lags+1:T-leads observations');
	else
		XX=[XX Ext];
	end
end

% OLS estimates: must check the non-singularity of (XX'XX): if its rank condition is below some
% tolerance value, the OLSE are computed by the QR factorisation as in REGRESS
[Q R]=qr(XX);
par=R\Q'*Y;
fit=XX*par;
res=Y-fit;
H=XX*(R\Q');
Y=center(Y);
ASS=diag(Y'*H*Y); % Adjusted sum of squares for each individual equation
TSS=diag(Y'*Y); % Total sum of squares for each individual equation
RSS=TSS-ASS;
%end;

con=sum((par(1:lags)));
ins=(par(lags+1));
