
function B=unweightvoting(A,class)

%%%-- Author: Shasha Mao, skymss@126.com(2011.11.21)---%%%
B=[];
m=size(A,1);
n=length(class);
%%% change the NaN to '0' %%%
for i=1:m
    asample=[];
    asample=A(i,:);
    indexNaN=[];
    indexNaN=find(isnan(asample));
    if (length(indexNaN~=0))
        asample(1,indexNaN)=0;
        A(i,:)=asample;
    end
end

%%% find the label that is not the class label of samples %%%
allclassA=unique(A);
none=length(allclassA);
number=zeros(none,1);
for i=1:none
    for j=1:n
        if (allclassA(i)~=class(j))
            number(i)=1;
        else
            number(i)=0;
        end
    end
end
indexnumber=find(number==1);
valuenoclass=[];
valuenoclass=allclassA(indexnumber,1);
%%% the number of no class %%%
numbernoclass=length(valuenoclass);
maxclass=zeros(m,n);   
for i=1:m       
    a=A(i,:);       
    for j=1:n       
        index=[];        
        index=find(a==class(j));          
        numberJ=length(index);          
        maxclass(i,j)=numberJ;     
    end
    vlauemax=[];      
    indexmax=[];   
    [valuemax(i,:),indexmax(i,:)]=max(maxclass(i,:));     
    B(i,1)=class(indexmax(i,:));  
end  
  



