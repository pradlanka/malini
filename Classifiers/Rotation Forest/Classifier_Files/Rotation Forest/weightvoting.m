
function [preEnClabel,sum]=weightvoting(preEachC,weightvalue,class)

%%%-- Author: Shasha Mao, skymss@126.com(2011.11.27)---%%%
preEnClabel=[];
[m,L]=size(preEachC);
numberclass=length(class);
sum=zeros(m,numberclass);
for i=1:m
    for j=1:L
        aprelabel=preEachC(i,j);
        aweightvalue=weightvalue(j,1);
        for k=1:numberclass            
            bvalue=HS(aprelabel,class(k));
            sum(i,k)=sum(i,k)+aweightvalue*bvalue;
        end
    end
end
[valuemax,indexmax]=max(sum,[],2);
preEnClabel=indexmax;

%%%%
function  value=HS(A,B)
if (A==B)
    value=1;
else
    value=0;
end