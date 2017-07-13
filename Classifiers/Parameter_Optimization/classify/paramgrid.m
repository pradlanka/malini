function A = paramgrid(Hyperparam)
% Modified from Allcomb to work for my code
NC = length(Hyperparam);

% check if we should flip the order
ii = NC:-1:1 ;

% check for empty inputs
if NC > 1
        % non-cell input, assuming all numerical values or strings
        % flip using ii if last column is changing fastest
        [A{ii}] = ndgrid(Hyperparam{ii}) ;
        % concatenate
        A = reshape(cat(NC+1,A{:}),[],NC) ;
elseif NC == 1
    A = Hyperparam{1}(:) ; % nothing to combine
end
