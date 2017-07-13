function [cum_perf_macro,cum_acc_macro,one_itr_acc] = cum_cal(cum_perf_macro, cum_acc_macro, performance, Confmattest)
v=isnan(performance (1,:,:)./(performance (1,:,:)+performance (4,:,:)));
for i=1:size(v,1)
    for j=1:size(v,2)
        for k=1:size(v,3)
           if v(i,j,k)
             vari(i,j,k)=0;
           else
            vari(i,j,k)= performance (1,j,k)./(performance (1,j,k)+performance (4,j,k));
           end
        end
    end
end
 cum_perf_macro(1,:,:) = cum_perf_macro(1,:,:)+ vari; % Precision
 cum_perf_macro(2,:,:) = cum_perf_macro(2,:,:)+ (performance (1,:,:)./(performance (1,:,:)+performance (2,:,:))); % Recall
 cum_perf_macro(3,:,:) = cum_perf_macro(3,:,:)+ (performance (3,:,:)./(performance (3,:,:)+performance (4,:,:))); % Specificity
 cum_acc_macro = cum_acc_macro+ (squeeze(sum(performance(1,:,:),2)./((sum(sum(Confmattest,1),2)))))'; % acc rate
 one_itr_acc= sum(performance(1,:,:),2)./(sum(sum(Confmattest,1),2)); % accuracy per iteration