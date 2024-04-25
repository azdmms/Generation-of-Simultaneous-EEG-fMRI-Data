function [finded,loc4]=choose_cluster(loc2,x,csize)
    %loc2=GridLoc(1:3:15003,:);
    difx=loc2-x;
    difx_norm=sum(difx.^2,2);
    sdif=sort(difx_norm);
    loc4=zeros(csize,3);
    finded=zeros(1,csize);
    for j=1:csize
       find1=find(difx_norm==sdif(j));
       loc4(j,:)=loc2(find1,:);
       finded(j)=find1;
    end
end
       
   
   