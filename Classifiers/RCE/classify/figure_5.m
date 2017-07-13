function figure_5(str,nsamples,ngenes,z_data,dist_i,s_names,filename,n_groups,unique_groups,groups)
methods = char('single','complete','average','weighted','centroid','median','ward');
dist = char('Euclidean','Correlation');
n = 7;
colors = char('red','green','blue','black','yellow');
marker = ['s' 'd' 'p' 'h' 'o'];
[n_methods mnmn] = size(methods);
maintitle = sprintf('%s #Top Significant Genes %d',str,ngenes );
y = pdist (z_data',dist(dist_i,:));
z = linkage(y,methods(n,:));

figure;
titl = sprintf('Hierarchical Cluster  Method:Inner squared distance(%s)','WARD' );
dendrogram(z,0,'colorthreshold','default','labels',s_names,'orientation','right');
dist_title = strcat('Pdist:',dist(dist_i,:),',',maintitle);
title(dist_title,'FontWeight','bold','FontSize',8,'interpreter','none');
[st,err] = sprintf('FileName:%s\n#Samples is %d  #Genes is %d',filename,nsamples,ngenes);

figure
dissimilarities = pdist(z_data',dist(dist_i,:));
[Y,stress] = mdscale(dissimilarities,2,'criterion','metricstress');
for j=1:n_groups
    hold on;
    grb  = unique_groups(j);
    loci = strcmp(groups,grb);
    plot(Y(loci,1),Y(loci,2),marker(j),'MarkerFaceColor',colors(j,:),'MarkerSize',10, 'MarkerEdgeColor','k');
end
text(Y(:,1),Y(:,2),s_names,'FontSize',8,'interpreter','none');
dist_title = strcat('Pdist:',dist(dist_i,:),',',maintitle);
title(dist_title,'FontWeight','bold');