function [accuracy,precision,recall,specificity]=acc_cal(cum_perf,FinalConfmat)
 accuracy = (squeeze(sum(cum_perf(1,:,:),2)./((sum(sum(FinalConfmat,1),2)))))';
 precision = squeeze((cum_perf(1,:,:)./(cum_perf(1,:,:)+cum_perf(4,:,:)))); % Precision
 recall = squeeze(cum_perf (1,:,:)./(cum_perf(1,:,:)+cum_perf(2,:,:))); % Recall
 specificity = squeeze(cum_perf(3,:,:)./(cum_perf(3,:,:)+ cum_perf(4,:,:))); % Specificity
 