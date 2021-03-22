function Edges=outToEdges(out)
Edges=[];
for i=1:size(out,1)
    a=out(i,:);
    a(a==0)=[];
    for j=2:length(a)
    Edges=[Edges;[a(1),a(j)]];
    end
end
        