
%Demonstrates the mcrvm algorithm on a toy example (Three class classification problem). Part of the code was
%derived from Tipping's matlab code.

%  @file mcrvm_example.m
%  				Author Arasanathan Thayananthan ( at315@cam.ac.uk)
%               (c) Copyright University of Cambridge
%  
%     This library is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 2 of the License, or (at your option) any later version.
%  
%     This library is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%  
%     You should have received a copy of the GNU Lesser General Public
%     License along with this library; if not, write to the Free Software
%     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function mcrvm_example

rand('state',1)

ALPHA_MAX		= 1e12;

synth=CreateData;

 
 N		= size(synth,1);
% P = number of classes in the data-1
 P=2; 
  
 


%Load the 2 class training data
X	= synth(:,1:2);
t	= zeros(N,P);

for p=1:P
    t(find(synth(:,3)==p),p)=1;
end
% Plot it
figure(1)
whitebg(1,'k')
clf
h_c1 = plot(X(synth(:,3)==1,1),X(synth(:,3)==1,2),'r.','MarkerSize',16);
hold on
h_c2 = plot(X(synth(:,3)==2,1),X(synth(:,3)==2,2),'y.','MarkerSize',16);
h_c3 = plot(X(synth(:,3)==3,1),X(synth(:,3)==3,2),'b.','MarkerSize',16);
box = 1.1*[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
axis(box)
drawnow

legend([h_c1 h_c2 h_c3],...
       'Class 1','Class 2','Class 3');
   title('3 class synthetic data','FontSize',14);


%set up the kernel matrix
kernel_	= '+gauss';
kernel_width		= 0.5;
PHI=sbl_kernelFunction(X,X,kernel_,kernel_width);
% the number of basis function (in this equal to N+1)
M=size(PHI,2); 


[probs, weights,alphas, used, kernels_]= mc_rvm(PHI,t,kernel_,3);



% showing results 

%
% Visualise the results
% 
gsteps		= 50;
range1		= box(1):(box(2)-box(1))/(gsteps-1):box(2);
range2		= box(3):(box(4)-box(3))/(gsteps-1):box(4);
[grid1 grid2]	= meshgrid(range1,range2);
Xgrid		= [grid1(:) grid2(:)];
% Compute RVM over a grid for visualisation purposes
% 

PHI2=cell(P,1);

for p=1:P
    PHI2{p}		= sbl_kernelFunction(Xgrid,X(used{p},:),kernels_{p},kernel_width);
end

y_grid=multinomial(PHI2,weights,size(Xgrid,1),P);

% Show decision boundary (p=0.5) and illustrate p=0.25 and 0.75
% 
[c,h05r]	= contour(range1,range2,reshape(y_grid(:,1),size(grid1)),[0.5 0.5],'-');
[c,h075r]	= ...
    contour(range1,range2,reshape(y_grid(:,1),size(grid1)),[0.25 0.75],'--');
set(h05r,'Color','r','LineWidth',3);
set(h075r,'Color',0.7*[1 0 0],'LineWidth',2);

[c,h05y]	= contour(range1,range2,reshape(y_grid(:,2),size(grid1)),[0.5 0.5],'-');
[c,h075y]	= ...
    contour(range1,range2,reshape(y_grid(:,2),size(grid1)),[0.25 0.75],'--');
set(h05y,'Color','y','LineWidth',3);
set(h075y,'Color',0.7*[1 1 0],'LineWidth',2);

[c,h05b]	= contour(range1,range2,reshape(1-y_grid(:,1)-y_grid(:,2),size(grid1)),[0.5 0.5],'-');
[c,h075b]	= ...
    contour(range1,range2,reshape(1-y_grid(:,1)-y_grid(:,2),size(grid1)),[0.25 0.75],'--');
set(h05b,'Color','b','LineWidth',3);
set(h075b,'Color',0.7*[0 0 1],'LineWidth',2);

for p=1:P
h_rv	= plot(X(used{p},1),X(used{p},2),'wo','LineWidth',2,'MarkerSize',10);
end
legend([h_c1 h_c2 h_c3 h05r h075r(1) h05y h075y(1) h05b h075b(1) h_rv],...
       'Class 1','Class 2','Class 3','Decision boundary','p=0.25/0.75','Decision boundary','p=0.25/0.75','Decision boundary','p=0.25/0.75','RVs')
   
hold off
title('RVM Classification of 3 class synthetic data','FontSize',14)



 function Y=multinomial(x,W,N,P)

Y=zeros(N,P);

sum=ones(N,1);
for p=1:P
    size(exp(x{p}*W{p}))
    sum=sum+ exp(x{p}*W{p});
end

for p=1:P
    Y(:,p)=exp(x{p}*W{p})./sum;
end