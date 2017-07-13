function figure_1(nc,no_classes,classname,accuracy,recall,specificity,precision,...
                 num_clusters,max_performance_each_cluster,avg_performance_each_cluster,...
                 TOP_RESAMPLE,upper_kfold,st,num)
             
    Daccuracy = accuracy;
    Drecall  = recall;
    Dspecificity = specificity;
    Dprecision = precision;

    classname=char(classname);
    figure;
    axes();
    set(gca,'XLim',[1 length(num_clusters)]);% This automatically sets the
    set(gca,'XTick',[1:1:length(num_clusters)]);  
    plot(gca,num,Daccuracy ,'-b','LineWidth',2,'Marker','h','MarkerFaceColor','y','MarkerSize',10); hold on
    plot(gca,num,Drecall ,'-r','LineWidth',2,'Marker','>','MarkerFaceColor','y','MarkerSize',10); hold on
    plot(gca,num,Dspecificity  ,'-g','LineWidth',2,'Marker','s','MarkerFaceColor','y','MarkerSize',10);
    plot(gca,num,Dprecision ,'-m','LineWidth',2,'Marker','^','MarkerFaceColor','y','MarkerSize',10);
    plot(gca,num,max_performance_each_cluster ,'-k','LineWidth',2,'Marker','x','MarkerFaceColor','y','MarkerSize',10);
        plot(gca,num,avg_performance_each_cluster ,'-m','LineWidth',2,'Marker','d','MarkerFaceColor','y','MarkerSize',10);
    set(gca,'XTick',[1:1:length(num_clusters)]);
    my_xticklabel_rotate(num,45,'',st,0.1,'FontSize',10);
    % legend('Average Performance','Max Performance','Combine Clusters Genes','%TP','%TN','LOCATION','BEST');
    legend('ACC of Combined Genes Clusters','%Precision','Specificity','Recall','Max','Avg','LOCATION','BEST');
    ylabel('Performance');
    xlabel('(Number of Clusers)/(Avg Number of Paths)');
    grid on;
    subtitle= sprintf('%s, Resample=%d,Kfold=%d Final Performance Across Clusters(micro)',classname,TOP_RESAMPLE,upper_kfold );
    title(subtitle);
hold on;