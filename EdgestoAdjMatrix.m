function Adj=EdgestoAdjMatrix(Edges)
ind=Edges(:,1:2);
ind=unique(ind);
ind=sort(ind);
for i=1:size(Edges,1)
    id1=Edges(i,1);
    id2=Edges(i,2);
    r1=find(ind==id1);
    r2=find(ind==id2);
    Adj(r1,r2)=1;
    Adj(r2,r1)=1;
end