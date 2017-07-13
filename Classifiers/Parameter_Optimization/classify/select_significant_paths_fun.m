function [mask] = select_significant_paths_fun(data_train_tmp,classname,groups_train)
%  Multivatiate N-way Anova (i.e. N-way hypothesis testing)  
    n=length(groups_train);
    grp=classname; 
    no_classes = length(classname);
    loc=zeros(no_classes,1); 
%     for  j=1:no_classes
%         for jj=1:n
%             (groups_train(jj)==grp(j))
%             loc(j)=loc(j)+(groups_train(jj)==grp(j));
%         end
%     end        
[groups_names,loc] = unique(groups_train);
loc =[loc; length(groups_train)+1];
    for j=1:no_classes
        d{j}=data_train_tmp(loc(j):loc(j+1)-1,:);
    end 
    if no_classes ==2
        [h,significance,ci] = ttest2(d{1},d{2},0.5); %filter genes based on the training set
    else
     significance = zeros(size(data_train_tmp,2),1);
     for ft=1:size(data_train_tmp,2)
         significance(ft)=anova1(data_train_tmp(:,ft),groups_train,'off');
     end
    end
    
    [pID,pN]= FDR(significance, 1); 
    mask = (significance<pID); %number of significant start paths by t-test after correcting for multiple comparisions
end
    
    
  function [pID,pN] = FDR(p,q)
% FORMAT [pID,pN] = FDR(p,q)
% 
% p   - vector of p-values
% q   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
%______________________________________________________________________________
% $Id: FDR.m,v 1.1 2009/10/20 09:04:30 nichols Exp $


%p = p(fininte(p));  % Toss NaN's
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID)));
pN = p(max(find(p<=I/V*q/cVN)));
    end
%end
