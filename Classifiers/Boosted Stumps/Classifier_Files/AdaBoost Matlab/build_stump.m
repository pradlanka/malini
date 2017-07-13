
% ----------------------------------------------
% build a stump from each component and return the best

function [stump] = build_stump(X,y,w)

d = size(X,2);
w = w/sum(w); % normalized the weights (if not already)

stump = cell(d,1); 
werr = zeros(d,1);
for i=1:d,
  stump{i} = build_onedim_stump(X(:,i),y,w);
  stump{i}.ind = i; 
  werr(i) = stump{i}.werr;
end;

[min_werr,ind] = min(werr);

stump = stump{ind(1)}; % return the best stump 

% ----------------------------------------------
% build a stump from a single input component

function [stump] = build_onedim_stump(x,y,w) 

[xsorted,I] = sort(x); % ascending 
Ir = I(end:-1:1); % descending

score_left  = cumsum(w(I).*y(I)); % left to right sums
score_right = cumsum(w(Ir).*y(Ir));  % right to left sums

% score the -1 -> 1 boundary between successive points 
score = -score_left(1:end-1) + score_right(end-1:-1:1); 

% find distinguishable points (possible boundary locations)

Idec = find( xsorted(1:end-1)<xsorted(2:end) );

% locate the boundary or give up

if (length(Idec)>0),
  [maxscore,ind] = max( abs(score(Idec)) ); % maximum weighted agreement
  ind = Idec(ind(1)); 

  stump.werr = 0.5-0.5*maxscore; % weighted error
  stump.x0   = (xsorted(ind)+xsorted(ind+1))/2; % threshold
  stump.s    = sign(score(ind)); % direction of -1 -> 1 change
else
  stump.werr = 0.5;
  stump.x0   = 0; 
  stump.s    = 1; 
end;
