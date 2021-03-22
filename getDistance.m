%% get distance between node nn(1) and nn(2) depending on whether they cut boundary or not
function len=getDistance(ps,Edges,E,Lx,Ly)
global rowID
rowID=getrowID(ps(:,3));
len=[];
E=sort(E);
i=find(Edges(:,1)==E(1) & Edges(:,2)==E(2));
if isempty(i)
    disp('not an edge in getDistance')
    pause
elseif length(i)==1
    r1=rowID(E(1));
    r2=rowID(E(2));
    if Edges(i,3)==0
         len = norm(ps(r1,1:2)-ps(r2,1:2));
    else
        psnew=ps(r2,1:2);
        % cut cross x
        if ps(r2,1)-ps(r1,1)>0.5*Lx
            psnew(1)=psnew(1)-Lx;
        elseif ps(r1,1)-ps(r2,1)>0.5*Lx        
            psnew(1)=psnew(1)+Lx;
        end  
        
        % cut cross y
        if ps(r2,2)-ps(r1,2)>0.5*Ly
            psnew(2)=psnew(2)-Ly;
        elseif ps(r1,2)-ps(r2,2)>0.5*Ly        
            psnew(2)=psnew(2)+Ly;
        end    
        len = norm(ps(r1,1:2)-psnew);  
    end
elseif length(i)>1
    disp('edge duplication')
    pause
end