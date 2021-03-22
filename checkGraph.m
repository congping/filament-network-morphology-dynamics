function flag=checkGraph(ps,d,out,Edges)
E=Edges(:,1:2);
[C,ic,~] = unique(E,'rows');
flag=true;
if size(E,1)~=size(C,1)
    flag=false;
    Edges=Edges(ic,:);
    disp('Edges duplicates');
    pause;
end

if sum(d)~=size(Edges,1)*2
    flag=false;
    disp('Edges row number is not consistent with degree')
    pause
end

ID1=sort(out(:,1));
ID2=sort(unique(Edges(:,1:2)));
if length(ID1)~=length(ID2)
    disp('ID in Edges not consistent with ID in out')
    flag=false;
    pause
else
    if ID1~=ID2
    disp('ID number in Edges not consistent with ID in out')
    flag=false;
    pause
    end
end




p=ps(:,1:2);
p1=unique(p,'rows');
if size(p)~=size(p1)
    flag=false;
    disp('ps position duplicates')
    pause
end



E=outToEdges(out);
E=sort(E,2);
E=sortrows(E);
E=unique(E,'rows');
E1=sort(Edges(:,1:2),2);
E1=sortrows(E1);
[a,b]=setdiff(E,E1,'rows');
if ~isempty(b)
    disp('Edges not consitent with out!')
    flag=false;
    pause
end



