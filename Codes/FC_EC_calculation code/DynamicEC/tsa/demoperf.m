% Demonstrates the much higher performance 
% of EARPYW.M compared to AR.M
% Both algoriths implemented the Yule-Walker approach
% for estimating the Autoregressive Parameters 

%	Version 2.43
%	last revision 19.06.1998
%	Copyright (c) 1997-1998 by Alois Schloegl
%	e-mail: a.schloegl@ieee.org	

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


if exist('OCTAVE_VERSION')==5
	 fprintf(2,'DEMOPERF cannot count FLOPS\n');
end;
% 	 fprintf(2,'DEMOPERF cannot be executed because FLOPS is not available under octave\n');
%        return;
%else
fprintf(1,'\nAR-test, calculates Autoregressive parameters with the Yule-Walker method\n');

clear eeg8s;load eeg8s.mat;
eeg8s=repmat(eeg8s,100,1);
MOP = 20;

if exist('ar')>=2
tic,flops(0); AR_ident=ar(eeg8s,MOP,'fb');RES=[flops toc];
which ar;
fprintf(1,'   needs \t%i FLOPS,%3.3fs for %i samples and a model order of %i.\n',RES,length(eeg8s),MOP);
else
        % AR.M is part of the Matlab IDENT-toolbox, which is not available 
        fprintf(1,'Sorry, AR.M not available \n');
end;        
tic,flops(0); AR_tsa=lattice(eeg8s',MOP);RES=[flops toc];
which lattice;
fprintf(1,'   needs \t%i FLOPS, %3.3fs for %i samples and a model order of %i.\n',RES,length(eeg8s),MOP);
    
fprintf(1,'\nHistogramm test\n');
tic,flops(0); [H,X]=hist(eeg8s,1060);RES=[flops toc];
which hist;
fprintf(1,'   needs \t%i FLOPS, %3.3fs for %i samples and a binwidth of 0.01 (1060 bins).\n',RES,length(eeg8s));
tmp=round(eeg8s/.01);
tic,flops(0); [R]=histo3(tmp);RES=[flops toc];
which histo;
fprintf(1,'   needs \t%i FLOPS, %3.3fs for %i samples and a binwidth of 0.01 (1060 bins).\n',RES,length(eeg8s));

fprintf(1,'\nUCP-Test: Test polynomials if all roots are inside the Unit Circle\n');
N=100;
x=rand(N,10);
y=zeros(N,11);

for k=1:N,y(k,:)=poly(x(k,:));end;

tic,flops(0); for k=1:N, b(k)=all(abs(roots(y(k,:)))<1); end; RES=[flops toc];
fprintf(1,'for k=1:N, b(k)=all(abs(roots(y(k,:)))<1); end;\n   needs \t%i FLOPS, %3.3fs for %i polynomials of order %i.\n',RES,N,MOP+1);

tic,flops(0); b=ucp(y);RES=[flops toc];
which ucp;
fprintf(1,'   needs \t%i FLOPS, %3.3fs for %i polynomials of order %i.\n',RES,N,MOP+1);

%end;
