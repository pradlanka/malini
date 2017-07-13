function figure_2(nc,classname,cum_acc_macro,cum_perf_macro,counter,num_paths,accuracy,recall,specificity,precision,tn_rate,...
                   num_clusters,max_performance_each_cluster,...
                   avg_performance_each_cluster,st_max,TOP_RESAMPLE,resample,upper_kfold,st,num)
figure
axes();
set(gca,'XLim',[1 length(num_clusters)]);
set(gca,'XTick',[1:1:length(num_clusters)]);
num = [1:1:length(num_clusters)];
plot(gca,num,max_performance_each_cluster ,'-k','LineWidth',2,'Marker','x','MarkerFaceColor','y','MarkerSize',10);
plot(gca,num,avg_performance_each_cluster ,'-m','LineWidth',2,'Marker','d','MarkerFaceColor','y','MarkerSize',10);
set(gca,'XTick',1:1:length(num_clusters));
my_xticklabel_rotate(num,45,'',st_max,0.1,'FontSize',10);
legend('Max','Avg','LOCATION','BEST');
ylabel('Performance');
xlabel('Avg Number of Paths');
classname=char(classname);
grid on;
subtitle = sprintf('%s, Resample=%d,Kfold=%d Final Performance Across Clusters(micro)',classname,TOP_RESAMPLE,upper_kfold );
title(subtitle);
res_out = [num_clusters;num_paths';accuracy;recall;specificity;precision;tn_rate]';
fout =  fopen('Finalaccurecy_micro.txt','w');
fprintf(fout,'Class=%s\n\n#Iteration\t=\t%d\n',classname,resample);
fprintf(fout,'#Clusters\t#paths\tAcc\tTP\tTN');
for resi=1:size(res_out,1)
    fprintf(fout,'\n%f\t%f\t%f\t%f\t%f',res_out(resi,1),res_out(resi,2),res_out(resi,3),res_out(resi,4),res_out(resi,5),res_out(resi,6),res_out(resi,7));
end
figure
precision = cum_perf_macro(1,:)./counter  ;
recall  = cum_perf_macro(2,:)./counter  ;
specificity = cum_perf_macro(3,:)./counter ;
accuracy = cum_acc_macro./counter ;
axes();
set(gca,'XLim',[1 length(num_clusters)]);% This automatically sets the
set(gca,'XTick',[1:1:length(num_clusters)]);
plot(gca,num,accuracy ,'-c','LineWidth',2,'Marker','d','MarkerFaceColor','y','MarkerSize',10); hold on
plot(gca,num,specificity ,'-b','LineWidth',2,'Marker','h','MarkerFaceColor','y','MarkerSize',10); hold on
plot(gca,num,recall ,'-r','LineWidth',2,'Marker','>','MarkerFaceColor','y','MarkerSize',10); hold on
plot(gca,num,precision ,'-g','LineWidth',2,'Marker','s','MarkerFaceColor','y','MarkerSize',10);
set(gca,'XTick',[1:1:length(num_clusters)]);
my_xticklabel_rotate(num,45,'',st,0.1,'FontSize',10);
legend('ACC of Combined paths Clusters','%Acc','%SP','%RE','PR','LOCATION','BEST');
% legend('accuracy','LOCATION','BEST');
ylabel('Performance');
xlabel('(Number of Clusers)/(Avg Number of Paths)');
grid on;
subtitle = sprintf('%s,Resample=%d,Kfold=%d Final Performance Across Clusters(macro)',classname,TOP_RESAMPLE,upper_kfold );
title(subtitle);