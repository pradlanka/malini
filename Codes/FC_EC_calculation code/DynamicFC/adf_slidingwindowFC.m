function [smooth_window_length, dynamicFC]=adf_slidingwindowFC(data,n4,n5)
% this method uses augmented Dickey-Fuller test (adf) plus temporal filtering to find sliding window length  
%n1 is #of timepoints, n2 is # of regions, n3 is # of subjects, data is
 %extracted time series from regions, and of dimension # of time by # of regions by
 %# of subjects
 %n4 is window length lower limit, n5 is window length upper limit, we
 %calculate window length with this range [n4,n5] using adf test. 
 %Output smooth_window_length is matrix of (n1-n5) by #of subjects,
 %dynamic_FC is dynamic Functional connectivity calculated, a matrix
 %(n1-n5) by # of regions by # of regions by %# of subjects

[n1,n2,n3]=size(data);
startpoint=n5;  %this is the start point to calculate dynamic FC, set to window length upper limit
 window_length=zeros(n1-startpoint,n3);smooth_window_length=zeros(n1-startpoint,n3);

 for k=1:n3  %loop over subjects
     k;
    vv=[];
    for i=(startpoint+1):n1  %loop over time
        v=n4;
        while v<startpoint;
            judge=zeros(n2,1);
        for j=1:n2
        result  =  adf(data( (i-v+1):  i, j,k) , -1, 1);   %adf test
             if result.adf >=result.crit(2); judge(j)=1;end   % result.crit(2) is for p=0.05 level
        end
        if judge ==0; 
            break; 
        else v=v+1;
        end;
        end
        vv=[vv; v]; 
    end
    window_length(:,k)=vv;
   
 end

 % FIR two pass smoothing of window length
 [B,A]=fir1(25,0.05); %we use FIR filter of order 25, 0.05 is cutoff frequency
  for k=1:n3,
      smooth_window_length(:,k)=round( filtfilt(B,A,window_length(:,k)));
  end
 
% dynamic FC    
dynamicFC=zeros(n1-startpoint,n2,n2,n3);  
for k=1:n3
      for i=(startpoint+1):n1   
            dynamicFC((i-startpoint),:,:, k) = corrcoef(  data((i-smooth_window_length(i-startpoint, k)+1): i ,:,k)); %using Pearson's correlation
      end
end
end
 
 


    
    
    
            
