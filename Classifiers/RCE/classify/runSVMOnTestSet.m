function [tp, fn, tn, fp, Confusion_mat,varargout]= runSVMOnTestSet(data,group_no,data_test,groups_test,test,classes,varargin)
% data: each arrow is one observation
% species : arrow(i) is the class member of abservation data(i)
% clsname : The name of the positive class or the case class

n_samples_tot= length(test);

samples_freq = zeros(1,n_samples_tot);% keep the frequency of the samples been tested
tp=0; tn=0; fp=0; fn=0;

% Change the group names to group indices for training and test dataset
cp = classperf(groups_test');
Decision_surf = classifier_train(data,group_no);
[classes_predict] = classifier_predict(Decision_surf, data_test);

% Updates the CP object with the current classification results
classperf(cp,classes_predict);
for j =1:length(test)
    if (test(j)==1)
        samples_freq(j)   = samples_freq(j)+1; %Keep the frequncy of the samples
    end
end

% Calculate performance Metrics
Con=confusionmat(groups_test,classes_predict,'order',classes);
L=length(unique(classes));

for KJ=1:L
    tp(KJ)=Con(KJ,KJ);
    fn(KJ)=sum(Con(KJ,:))-Con(KJ,KJ);
    fp(KJ)=sum(Con(:,KJ))-Con(KJ,KJ);
    tn(KJ)=sum(sum(Con))-sum(Con(KJ,:))-sum(Con(:,KJ))+Con(KJ,KJ);
end
Confusion_mat=Con;
if nargin > 6
Decision_Surface_Flag =varargin{1};
if Decision_Surface_Flag == 1
    varargout{1}= Decision_surf;
end 
end
