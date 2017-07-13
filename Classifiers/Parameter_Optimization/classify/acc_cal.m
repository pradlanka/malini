function [precision,recall,specificity]=acc_cal(FinalConfmat)
precision =zeros(1,size(FinalConfmat,1)); recall = zeros(1,size(FinalConfmat,1)); specificity = zeros(1,size(FinalConfmat,1));
for i=1:size(FinalConfmat,1)
 precision(i) = FinalConfmat(i,i)./ sum(FinalConfmat(:,i));% Precision
 recall(i) = FinalConfmat(i,i)/ sum(FinalConfmat(i,:)); % Recall
 tn = sum(sum(FinalConfmat)) - sum(FinalConfmat(i,:))-sum(FinalConfmat(i,:))+ FinalConfmat(i,i); 
 fp = (sum(FinalConfmat(:,i))- FinalConfmat(i,i));
 specificity(i) = tn/(tn+fp); % Specificity
end