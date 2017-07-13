
function conn_thresh=threshold(conmat,conn_sur,num_std,option)

% conmat and conn_sur from main_cpgc_func

% num_std is the number of standard deviations away from the mean. Use 2 or
% 3. You can increase the threshold by increasing this number

% option=0 for regular connectivity matrix and 1 for difference
% connectivity matrix

if option==0 %regular
    
for i=1:size(conmat,1)
    for j=1:size(conmat,1)
        if i~=j
            if conmat(i,j)>=mean(conn_sur(i,j,:))+num_std*std(conn_sur(i,j,:))
                conn_thresh(i,j)=conmat(i,j);
            else
                conn_thresh(i,j)=0;
            end
        else
            conn_thresh(i,j)=0;    
        end
    end
end

end

if option==1 % difference
    
for i=1:size(conmat,1)
    for j=1:size(conmat,1)
        if i~=j
            if conmat(i,j)>0
                
            if conmat(i,j)>=mean(conn_sur(i,j,:))+num_std*std(conn_sur(i,j,:))
                conn_thresh(i,j)=conmat(i,j);
            else
                conn_thresh(i,j)=0;
            end
            
            end
            
            if conmat(i,j)<0
                
            if conmat(i,j)<=mean(conn_sur(i,j,:))-num_std*std(conn_sur(i,j,:))
                conn_thresh(i,j)=conmat(i,j);
            else
                conn_thresh(i,j)=0;
            end
            
            end
            
        else
            conn_thresh(i,j)=0;    
        end
    end
end

end