function [R,tix]=histo3(Y)
% HISTO3 calculates histogram and performs data compression
% 
% R = HISTO3(Y)
% 	R is a struct with th fields 
%       R.X  are the bin-values 
%       R.H  is the frequency of occurence of value X 
%  	R.N  are the number of valid (not NaN) samples 
%
% Data compression can be performed in this way
%   	[R,tix] = histo3(Y) 
%      		is the compression step
%
%	R.tix provides a compressed data representation. 
%	R.compressionratio estimates the compression ratio
%
% 	R.X(tix) and R.X(R.tix) 
%		reconstruct the orginal signal (decompression) 
%
% The effort (in memory and speed) for compression is O(n*log(n)).
% The effort (in memory and speed) for decompression is O(n) only. 
%
% see also: HISTO, HISTO2, HISTO3, HISTO4
%
% REFERENCE(S):
%  C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).


%	$Id: histo3.m 5090 2008-06-05 08:12:04Z schloegl $
%	Copyright (C) 1996-2002,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/
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


[yr,yc]=size(Y);
if yr==1,
        % Makes sure there is a second row
        % Sort does not support the DIM-argument, therefore,
        % this function would not work correctly with this software
        % Once this is fixed, this part can be removed. 
        Y = [Y; NaN+ones(size(Y))];  
end;

% identify all possible X's and overall Histogram
[sY ,idx] = sort(Y(:));
[tmp,idx] = sort(idx);        % generate inverse index

ix  = diff(sY,1)>0;
tmp = [find(ix); sum(~isnan(sY))];
H   = diff([0; tmp]);

R.datatype = 'HISTOGRAM';
R.X = sY(tmp);
R.N = sum(~isnan(Y),1);


% generate inverse index
if nargout>1,
        tix = cumsum([1;ix]);	% rank 
        tix = reshape(tix(idx),yr,yc);		% inverse sort rank
        cc  = 1;
        tmp = sum(ix)+1;
	if exist('OCTAVE_VERSION')>=5,
		; % NOP; no support for integer datatyp 
        elseif tmp <= 2^8;
                tix = uint8(tix);
                cc = 8/1;
        elseif tmp <= 2^16;
                tix = uint16(tix);
                cc = 8/2;
        elseif tmp <= 2^32;
                tix = uint32(tix);
                cc = 8/4;
        end;
        R.compressionratio = (prod(size(R.X)) + (yr*yc)/cc) / (yr*yc);
	R.tix = tix;        
end;


% if yc==1, we are all set; else 
if yc>1,	% a few more steps are necessary
        H0 = H; %overall histogram	
        % allocate memory
        H = zeros(size(R.X,1),yc);
        
        % scan each channel
        for k = 1:yc,
		sY = sort(Y(:,k));
		ix = find(diff(sY,1)>0);
                if size(ix,1)>0,
                        tmp = [ix; R.N(k)];
                else
                        tmp = R.N(k);
                end;

                t = 0;
                j = 1;
                for x = tmp',
                        acc = sY(x);
                        while R.X(j)~=acc, j=j+1; end;
                        %j = find(sY(x)==R.X);   % identify position on X 
                        H(j,k) = H(j,k) + (x-t);  % add diff(tmp)
                        t = x;
                end;
        end;
        
        if any(H0~=sum(H,2)),  %%% CHECK 
                fprintf(2,'ERROR HISTO\n');
        end;	
end;

R.H = H;


 