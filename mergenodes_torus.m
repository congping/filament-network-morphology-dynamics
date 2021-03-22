function [pn,dn,outn,Ntn,Edgesn,flag]=mergenodes_torus(ps,d,out,Nt,Lx,Ly,Edges)

%merge nodes if their distance is less than eta  
% out: [i, j1,j2,...], ps - positions, d- degree for the same node in the
% same row
global eta Npp rowID
flag=false;

%Ntn=Nt;
pn=ps;
dn=d;
outn=out;
Edgesn=Edges;

Edges=sortrows(Edges,1);
D=Edges(:,4);%D=pdist(ps);

ds=sort(find(D<eta),'descend');
rdelete=[];
for i=1:length(ds)
    rowID=getrowID(out);
    nn1=Edges(ds(i),1);
    nn2=Edges(ds(i),2);  
    r1=rowID(nn1);r2=rowID(nn2);
    rdelete=[rdelete,r2];
    outn(find(outn==nn2))=nn1;

    %%% care of duplicated edges and check cut boundary
    % deleting edges link to nn2, and creat new edges link to nn1
    for k=1:dn(r2)
        id1=outn(r2,k+1);
        if id1~=nn1 && id1~=0
            [~,f1,~]=checkCutBoundary(Edgesn,id1,nn1);
            if ~f1 % if not an edge then add one Edge        
                    r3=rowID(id1);    
                    [d3, flag3]=getDistanceCut(pn,r1,r3,Lx,Ly);
                    Edgesn(end+1,:)=[sort([id1,nn1]),flag3,d3]; 
            else
                    disp('not an edge in merging');
            end
        end
    end
    ind=find(Edgesn(:,1)==nn2);
    Edgesn(ind,:)=[];
    ind=find(Edgesn(:,2)==nn2);
    Edgesn(ind,:)=[];
    
    
    
% move neighbour of nn2 to nn1
neigh2=outn(r2,2:end);
neigh1=outn(r1,2:end);
neigh=setdiff(union(neigh1,neigh2),[nn1,0]);% remove 0
dn(r1)=length(neigh);
outn(r1,2:end)=0;
outn(r1,2:2+dn(r1)-1)=neigh;
outn(r2,:)=0;% set to zero in order to avoid duplicate of node ID in the first collum

%replacing ID nn2 by node nn1

% repack (remove duplicate IDs) neighors of nn1, nn2
for j=1:dn(r1)
    rj=find(outn(:,1)==neigh(j));
    nei=setdiff(outn(rj,2:end),[neigh(j),0]);
    outn(rj,2:end)=0;
    outn(rj,2:2+length(nei)-1)=nei;
    dn(rj)=length(nei);
end

if ~ismember(nn1,outn(:,1))
    display('error: node not in out');
   break
end
%display(strcat('merge node ',num2str(nn2),'to node ',num2str(nn1)));

end
% merge process
%rdelete
pn(rdelete,:)=[];
outn(rdelete,:)=[];
dn(rdelete)=[];
Ntn=Nt-length(unique(rdelete));

if ~isempty(ds)
    flag=true;
    %checkGraph(pn,dn,outn,Edgesn);
end