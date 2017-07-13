%  @file mc_rvm_example.m
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


%This m file implements the bottom-up  relevance vector machine for
%multi-class RVM classification , details can be found in the following
%documents

%   Multivariate Relevance Vector Machines for Tracking 
%   Arasanathan Thayananthan et al. (University of Cambridge)
    
 %   in Proc. 9th European Conference on Computer Vision 2006.

 %Detailed derivation can be found in 
 %Arasanathan Thayananthan. Template-based pose estimation and tracking of 3 hand motion
 %PhD thesis, University of Cambridge, UK 2005.

% Relevance Vector Machine based Mixture of Experts 
%A. Thayananthan, Technical Report, Department of Engineering, University of Cambridge, September 2005.


function [probs, weights,alphas, used, chosen_kernels_]= mc_rvm(PHI,tdata,kernel_,maxIts)

% dimensions
% N number of training data
% P number of classes - 1
% M numer of basis functions

% Input
% Matrix PHI   - N x M - kernel matrix
% Matrix tdata - N x P - classfication matrix e.g.  for a four class
                           % problem it will look like [0 1 0;0 0 1 ;0 0 1;1 0 0......]
%kernel_ kernel type- to handle bias.
 % maxIts - maximum number of iterations                         


%Output

% Matrix  probs - N x P -  final classiffication probability for each training data for each
% 
% Cell of vector  weights - Cell(P,1)- weights for chosen basis functions for each class
% Cell of vector  alphas -  Cell(P,1)- alpha values for each basis function for each class
% chosen kernels - conatins information of whether bias is chosen for each class






ALPHA_MAX		= 1e12;
P=size(tdata,2);
[N M]=size(PHI);

chosen_kernels_=cell(P,1);


alphas=cell(P,1);
weights=cell(P,1);
for p=1:P
    alphas{p}=zeros(M,1);
    alphas{p}(:)=ALPHA_MAX;
    alphas{p}(1)=1.0;
    weights{p}=zeros(1,1);
end




[weights used  alphas] = UpdateParams(PHI,tdata,weights,alphas,maxIts);

PHI2=cell(P,1);

for p=1:P
    chosen_kernels_{p}=kernel_;
    
    PHI2{p}=PHI(:,used{p});
    
    if  chosen_kernels_{p}(1)=='+'	
  % Take account of bias if originally used ...
        used{p}	= used{p} - 1;
        if used{p}(1)~=0
        % ... and if pruned ...
            chosen_kernels_{p}(1)	= [];
        else
            used{p}(1)	= [];
        end
   
    end
    
    
end
 
probs=multinomial(PHI2,weights,N,P);


function Y=multinomial(x,W,N,P)

Y=zeros(N,P);

sum=ones(N,1);
for p=1:P
    sum=sum+ exp(x{p}*W{p});
end

for p=1:P
    Y(:,p)=exp(x{p}*W{p})./sum;
end


function [ weights, used, alpha] =UpdateParams(PHI,T,weights,alpha,maxIts)

ALPHA_MAX		= 1e12;
[N,M]	= size(PHI);
P=size(T,2);
nz=zeros(P,1);
nonZero=cell(P,1);
alpha_nz=cell(P,1);
w=cell(P,1);
w_nz=cell(P,1);
PHI_nz=cell(P,1);
for p=1:P
    nonZero{p}	= (alpha{p}<ALPHA_MAX);
    w{p}=zeros(M,1);
    w{p}(nonZero{p})=weights{p};
end



for i=1:maxIts

  for p=1:P
    nonZero{p}	= (alpha{p}<ALPHA_MAX);
    alpha_nz{p}	= alpha{p}(nonZero{p});
    w{p}(~nonZero{p})	= 0;
    nz(p)=size(alpha_nz{p},1);
    w_nz{p}=w{p}(nonZero{p});
    PHI_nz{p}	= PHI(:,nonZero{p});
    
  end
  


  [w_nz  SIGMA_nz betaED] = EstimateWeights(PHI_nz,T,alpha_nz,w_nz,25);

  
    for p=1:P
        w{p}(nonZero{p})=w_nz{p};
        w_nz{p};
    end
    logBeta	= 0;
    Y=multinomial(PHI_nz,w_nz,N,P);
    
    B=zeros(N,N,P);
    t_hat=zeros(N,P);
    for p=1:P
        B(:,:,p)=diag(Y(:,p).*(1-Y(:,p)));
        t_hat(:,p)=PHI_nz{p}*w_nz{p}+inv(B(:,:,p))*(T(:,p)-Y(:,p));
    end
   
    
    max_change=-inf;
    max_k=-1;
    max_alpha=-inf;
    
   
    
   for p=1:P
   
    for k=1:M
         phi=PHI(:,k);
        for j=1:P
            S(j)=phi'*B(:,:,j)*phi-phi'*B(:,:,j)*PHI_nz{p}*SIGMA_nz{p}*PHI_nz{p}'*B(:,:,j)*phi;
            Q(j)=phi'*B(:,:,j)*t_hat(:,j)- phi'*B(:,:,j)*PHI_nz{p}*SIGMA_nz{p}*PHI_nz{p}'*B(:,:,j)*t_hat(:,j);
        end
        
        if(alpha{p}(k)<ALPHA_MAX)
            s=alpha{p}(k)*S./(alpha{p}(k)-S);
            q=alpha{p}(k)*Q./(alpha{p}(k)-S); 
        else 
            s=S;
            q=Q;
        end
        
        [new_alpha,l_inc]=SolveForAlpha(s',q',S',Q',alpha{p}(k));
        l_inc_vec(k)=l_inc;
        alpha_vec(k)=new_alpha;
        
                                                                                         
        if(max_change<l_inc)
            max_k=k;
            max_s=s;
            max_q=q;
            max_alpha=new_alpha;
            max_change=l_inc;
        end
    end
   
   if(max_alpha<ALPHA_MAX)
        alpha{p}(max_k)=max_alpha;
    else 
        break;
    end
 
  end
end

 for p=1:P
    nonZero{p}	= (alpha{p}<ALPHA_MAX);
    alpha_nz{p}	= alpha{p}(nonZero{p});
    w{p}(~nonZero{p})	= 0;
    nz(p)=size(alpha_nz{p},1);
    w_nz{p}=w{p}(nonZero{p});
    PHI_nz{p}	= PHI(:,nonZero{p});
    
 end

[w_nz  SIGMA_nz betaED] = EstimateWeights(PHI_nz,T,alpha_nz,w_nz,25);
% Tidy up return values
weights=cell(P,1);
used=cell(P,1);

for p=1:P
    weights{p}=w_nz{p};
    used{p}=find(nonZero{p});
end




function [w, SIGMA_nz, betaED] = EstimateWeights(PHI,T,alpha,w,its)


STOP_CRITERION	= 1e-6;
LAMBDA_MIN	= 2^(-8);

[P dd]	= size(PHI);
d=zeros(P,1);
for p=1:P
    [N d(p)]=size(PHI{p});
end
[M dd]	= size(w{1});
g=cell(P,1);

U=cell(P,1);
SIGMA_nz=cell(P,1);

w_new=cell(P,1);
delta_w=cell(P,1);

errs	= zeros(its,1);

Y	= multinomial(PHI,w,N,P);

data_term=0;
regulariser=0;

ss=min(size(find(Y==0)));
if(ss==0)

    for p=1:P
        data_term =data_term- sum(T(:,p).*log(Y(:,p)))/N;
        regulariser=regulariser+(alpha{p}'*(w{p}.^2))/(2*N);
    end
else 
    data_term=inf;
end
ss=min(size(find(sum(Y,2)==1)));
if(ss==0)
    data_term=data_term-sum((1-sum(T,2)).*log(1-sum(Y,2)))/N;
else 
    data_term=inf;
end

for p=1:P
       regulariser=regulariser+(alpha{p}'*(w{p}.^2))/(2*N);
end
    

err_new		=  data_term + regulariser;

break_condition=0;



for i=1:its
    
    for p=1:P
       vary	= Y(:,p).*(1-Y(:,p));
       
       PHIV	= PHI{p} .* (vary * ones(1,d(p)));
       e	= (T(:,p)-Y(:,p));
        A=diag(alpha{p});
        
  
        g{p}		= PHI{p}'*e - alpha{p}.*w{p};
        
        Hessian	= (PHIV'*PHI{p} + A);
        SIGMA_nz{p}=inv(Hessian);
        U{p}=chol(Hessian);
      
      if i==1
        condHess	= rcond(Hessian);
        if condHess<eps
%           fprintf(2,'(postMode) warning: ill-conditioned Hessian (%g)\n', ...
% 		      condHess);
%           fprintf(2,'(postMode) returning immediately for alpha-update\n');
          return
        end
      end
      
   
    end  
      
  
   errs(i)	= err_new;
   
%     
%    fprintf('PostMode Cycle: %2d\t error: %.6f\n', i, errs(i));
%        
   
    if( i>=2)
        break_condition=1;
        for p=1:P
            if (norm(g{p})/M>=STOP_CRITERION)
                break_condition=0;
            end
        end
    end
          
    
    if (break_condition==1)
        errs	= errs(1:i);
      
%           fprintf(['(postMode) converged (<%g) after %d iterations, ' ...
% 		    'gradient = %g\n'], STOP_CRITERION,i,norm(g{1})/M);
      
        break
   end
 
  for p=1:P
    delta_w{p}	= U{p} \ (U{p}' \ g{p});
    
  end
  lambda	= 1;
  while lambda>LAMBDA_MIN
      for p=1:P
        w_new{p}	= w{p} + lambda*delta_w{p};
      end
    Y=multinomial(PHI,w_new,N,P);
    data_term=0;
    regulariser=0;
   
    
    ss=min(size(find(Y==0)));
    if(ss==0)

        for p=1:P
            data_term =data_term- sum(T(:,p).*log(Y(:,p)))/N;
            regulariser=regulariser+(alpha{p}'*(w{p}.^2))/(2*N);
        end
    else 
        data_term=inf;
    end
    ss=min(size(find(sum(Y,2)==1)));
    if(ss==0)
        data_term=data_term-sum((1-sum(T,2)).*log(1-sum(Y,2)))/N;
    else 
        data_term=inf;
    end

    for p=1:P
       regulariser=regulariser+(alpha{p}'*(w{p}.^2))/(2*N);
    end
    
    
    err_new		=  data_term + regulariser;
     
    if err_new>errs(i)
      lambda	= lambda/2;
     
% 	    fprintf(['(postMode) error increase! Backing off ... (' ...
% 	    ' %.3f)\n'], lambda);
      
    else
      w		= w_new;
      lambda	= 0;
    end
 end
  if lambda
%     
%       fprintf(['(postMode) stopping due to back-off limit,' ...
% 			 'gradient = %g\n'],sum(abs(g{1})));
   
    break;
  end
end


%betaED	= e'*(e.*vary);
betaED=0;

%%
%%




function [new_alpha,l_inc]=SolveForAlpha(s,q,S,Q,old_alpha)
ALPHA_MAX		= 1e12;

[n m]=size(s);

index=[1:n]';
C=zeros(n,1);
CC=zeros(2*n-1,1);
for i=1:n
   SS=-s(find(index~=i));
   P=poly(SS)';
   PP=conv(P,P);
   CC=CC+q(i)*q(i)*PP;
   C=C+P;
end
P=poly(-s);
PP=n*conv(P,P)';
CC0=conv(P,C)';
CC=[0
    CC];

CCC=[CC+CC0
    0];

r=roots(PP-CCC);



r=r(find(abs(imag(r))==0));
r=r(find(r>0));
r=[r
    ALPHA_MAX];
[nn mm]=size(r);

L=zeros(nn,1);
for i=1:nn
    new_alpha=r(i);
    if((old_alpha>=ALPHA_MAX)& (r(i)<ALPHA_MAX))
        tt(i)=1;
        for j=1:n
            L(i)=L(i)+ log(new_alpha/(new_alpha+s(j)))+((q(j)*q(j)./(new_alpha+s(j))));
        end
    elseif((old_alpha<ALPHA_MAX)&(new_alpha>=ALPHA_MAX))
        tt(i)=2;
        for j=1:n
            L(i)=L(i)+ (Q(j)*Q(j)/(S(j)-old_alpha))-log(1- (S(j)/old_alpha));
        end
    elseif((old_alpha<ALPHA_MAX)&(new_alpha<ALPHA_MAX))
        tt(i)=3;
        for j=1:n
            ad=(1/new_alpha)-(1/old_alpha);
            L(i)=L(i)+ (Q(j)*Q(j)/(S(j)+(1/ad))) -log(1+S(j)*ad);
        end
    else
        tt(i)=4;
        L(i)=0;
    end
end
[m ind]=max(L);
new_alpha=r(ind);
l_inc=L(ind);