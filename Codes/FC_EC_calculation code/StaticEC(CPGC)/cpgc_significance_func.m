
function [conmat,p_thresh,p_val]=cpgc_significance_func(data,ord)

% data is a matrix which has all the time series stacked into it as columns
% ord is model order

for i=1:size(data,2)
        for j=1:size(data,2)
            if i~=j
                conmat(i,j)=CPGC_cap_func(data(:,i),data(:,j),ord);
            else
                conmat(i,j)=0;
            end
        end
end
    


for k=1:1000 % reduce value of k if it takes too much time
    surrogate=gen_surrogate(data,k);
    for i=1:size(data,2)
        for j=i:size(data,2)
            if i~=j
                conn_sur(i,j,k)=CPGC_cap_func(surrogate(20:size(data,1)-20,i),surrogate(20:size(data,1)-20,j),ord);
                conn_sur(j,i,k)=CPGC_cap_func(surrogate(20:size(data,1)-20,j),surrogate(20:size(data,1)-20,i),ord);
            else
                conn_sur(i,j,k)=0;
            end
        end
    end
    clear surrogate; 
end

for i=1:size(data,2)
    for j=1:size(data,2)
        if i~=j
            start=min(conn_sur(i,j,:));
            stop=max(conn_sur(i,j,:));
            vec=start:0.0001:stop;
            [n,x]=hist(squeeze(conn_sur(i,j,:)),vec); %%%% behnaz: I added the squeeze to remove error
                
            for m=1:length(x)
                if x(m)>=conmat(i,j)
                    p(i,j)=trapz(n(m:length(x)))/trapz(n);break
                end
            end
          
        end
    end
end

% If it returns very few significant paths, there is a way to do a less
% conservative test by using the following loop instead of the loop above

%for i=1:size(data,2)
%    for j=1:size(data,2)
%        if i~=j
%            [h,p(i,j)]=ttest(squeeze(conn_sur(i,j,:),conmat(i,j)));
%        end
%    end
%end
p_val = p;

for i=1:size(data,2)
    for j=1:size(data,2)
        if i~=j
            if p(i,j)<0.05
                p_thresh(i,j)=1;
            else
                p_thresh(i,j)=0;
            end
        end
    end
end





