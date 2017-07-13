function savevars(name,varargin)
numvars = (nargin-1)/2;
 for i = 1:numvars
     fname{i} = varargin{(2*i)-1};
     eval([fname{i},' = varargin{2*i};']);  
 end
save('-mat',name,fname{1});
for i = 2:numvars
    save('-mat',name,fname{i},'-append');
end
end