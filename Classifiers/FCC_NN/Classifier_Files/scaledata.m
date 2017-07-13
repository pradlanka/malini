function [dataout,varargout] = scaledata(datain,varargin)
if nargin == 3 
    datamin = varargin{1};
    scale = varargin{2};
    diff = datain - datamin;
else
datamin = min(datain(:));
diff = datain - datamin;
scale = 2*range(diff(:));
end
dataout = diff/scale;
if nargin == 1 && nargout ==3
    varargout{1} = datamin;
    varargout{2} = scale;
end
