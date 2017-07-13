
% this method uses augmented Dickey-Fuller test (adf) plus temporal filtering to find sliding window length  
%n1 is #of timepoints, n2 is # of regions, n3 is # of subjects, data is
 %extracted time series from regions, and of dimension # of time by # of regions by
 %# of subjects
 %n4 is window length lower limit, n5 is window length upper limit, we
 %calculate window length with this range [n4,n5] using adf test. 
 %Output smooth_window_length is matrix of (n1-n5) by #of subjects,
 %dynamic_FC is dynamic Functional connectivity calculated, a matrix
 %(n1-n5) by # of regions by # of regions by %# of subjects
 %load PTSD_All.mat
[n1,n2,n3]=size(PTSD);
startpoint=20;  %this is the start point to calculate dynamic FC, set to window length upper limit
 window_length=20*ones(n1-startpoint,n3); smooth_window_length=20*ones(n1-startpoint,n3);   
 
% dynamic FC    
dynamicFC =zeros(n1-startpoint,n2,n2,n3);  
for k=1:n3
      for i=(startpoint+1):n1   
            dynamicFC((i-startpoint),:,:, k) = corrcoef(PTSD((i-smooth_window_length(i-startpoint, k)+1): i ,:,k)); %using Pearson's correlation
      end
end
meanDFC=squeeze(mean(dynamicFC,1));
save('meanDFC_PTSD','meanDFC');
size(size(size(isnan(meanDFC),1),2),3)

            
