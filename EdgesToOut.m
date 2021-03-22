function [out,d]=EdgesToOut(Nt,Edges)
d=[];
out=[];
edg=[Edges(:,1:2);Edges(:,2:-1:1)];
for i=1:Nt %node id
    ind=find(edg(:,1)==i);
    %need to change to rowID
    d(i)=length(ind);
    out(i,1:length(ind)+1)=[i,edg(ind,2)'];
end