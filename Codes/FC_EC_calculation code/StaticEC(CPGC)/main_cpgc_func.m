function [conmat]=main_cpgc_func(data,order)
% data is extracted time series from regions, size(data)=M x N , where M=number of time points, N is number of regions
% The time series must be detrended and zscored (or normalized)

%order is the model order. For low pass filtered resting state data, I
%would try 1 to 5. 

%Calculation of CPGC (Static EC) into conmat, which is N x N.
for k=1:size(data,3)
for i=1:size(data,2)
        for j=1:size(data,2)
            if i~=j
                conmat(i,j,k)=CPGC_func(data(:,i),data(:,j),order);
            else
                conmat(i,j,k)=0;
            end
        end
end
end
%Calculation of null distribution of static EC using surrogate
%data,conn_sur, generate 1000 surrogate data of dimension N x N.
% for k=1:1000
%     surrogate=gen_surrogate(data,k);
%     for i=1:size(data,2)
%         for j=i:size(data,2)
%             if i~=j
%                 conn_sur(i,j,k)=CPGC_func(surrogate(20:size(data,1)-20,i),surrogate(20:size(data,1)-20,j),1);
%                 conn_sur(j,i,k)=CPGC_func(surrogate(20:size(data,1)-20,j),surrogate(20:size(data,1)-20,i),1);
%             else
%                 conn_sur(i,j,k)=0;
%             end
%         end
%     end
%     clear surrogate; 
% end



