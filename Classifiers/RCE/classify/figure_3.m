function figure_3(groups,groups_train,num_genes,samples_scores,total_freq,...
                    uniq_groups,n_samples,s_names,TOP_RESAMPLE,upper_kfold,...
                    tp_rate,tn_rate,accuracy);
fgrp = groups(1);
i=2:n_samples
sep=sum(ismember(groups(i),fgrp))+1;
l_indx = length(num_genes);
num_plots = 6;
for i=0:num_plots
    figure
    if (i == num_plots) 
        sum_scores = sum(samples_scores)'./(total_freq);
    else
        sum_scores = samples_scores(l_indx-i,:);
    end
    nelements=length(sum_scores);
    y1=zeros(length(sum_scores),1)';
    y2=y1;
    y1(1:sep)=sum_scores(1:sep);
    y2(sep+1:nelements)=sum_scores(sep+1:nelements);
    h = bar(y1,'r');
    hold on;
    bar(y2,'g');
    legend(uniq_groups,'Location','BEST','interpreter','none');
    my_xticklabel_rotate([1:1:n_samples],90,'',s_names,0.3);
    if (i==num_plots) 
        subtitle= sprintf('Samples Behavior all Over the Iteration','FontSize',10 );
    else
        subtitle=sprintf('Resample=%d,Kfold=%d,#Avg Genes=%5.2f Tp=%2.1f%% Tn=%2.1f%% Acc=%2.1f%%',...
        TOP_RESAMPLE,upper_kfold,num_genes(l_indx-i),tp_rate(l_indx-i)*100,...
        tn_rate(l_indx-i)*100,accuracy(l_indx-i)*100 );
    end
    title(subtitle);
end