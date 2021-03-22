function rowID=getrowID(out)
rowID=zeros(max(out(:,1)),1);
%out
for i=1:size(out,1)
	%i
	%size(out,1)
	%find(out(:,1)==out(i,1))
    rowID(out(i,1))=find(out(:,1)==out(i,1));
end
