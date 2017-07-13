function R=histo2(Y)
% HISTO2 calculates histogram of each column
% R=HISTO2(Y)
%
% R=HISTO(...)            
% 	R is a struct with th fields 
%       R.X  are the bin-values 
%       R.H  is the frequency of occurence of value X 
%  	R.N  are the number of valid (not NaN) samples 
%
% more histogram-based results can be obtained by HIST2RES2  
%
% see also: HISTO, HISTO2, HISTO3, HISTO4
%
% REFERENCE(S):
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).

%	$Id: histo2.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the TSA-toolbox 
%	http://hci.tugraz.at/~schloegl/matlab/tsa/
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


[yr,yc] = size(Y);

if yr==1,
        % Makes sure there is a second row
        % Sort does not support the DIM-argument, therefore,
        % this function would not work correctly with this software
        % Once this is fixed, this part can be removed. 
        Y = [Y; NaN+ones(size(Y))];  
end;

sY  = sort(Y);
N   = sum(~isnan(Y),1);
[ix,iy] = find(diff(sY,1)>0);
nn0 = 0;
warning('off');

for k=1:yc,
	tmp    = [ix(iy==k); N(k)];
        nn1    = length(tmp);
        
        H(1:nn1,k) = diff([0; tmp]);
       	X(1:nn1,k) = sY(tmp,k);

        if nn1<nn0;
                H(1+nn1:nn0,k) = NaN;
                X(1+nn1:nn0,k) = NaN;
        elseif nn1>nn0,
                H(1+nn0:nn1,1:k-1) = repmat(NaN,nn1-nn0,k-1);
                X(1+nn0:nn1,1:k-1) = repmat(NaN,nn1-nn0,k-1);
                nn0 = nn1;
        end;
end;

R.datatype = 'HISTOGRAM';
R.H = H;
R.X = X;
R.N = sumskipnan(R.H,1);


 