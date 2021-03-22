function [dn,outn,Edgesn,flag]=checkintersection_torus(ps,d,out,Lx,Ly,Edges,XY1,XY2)
% this is checking potential edge intersection, if so, merge intersection
% point to neares node. So node number does not change, only the topolgy
% changes, out degree would increase for some nodes
global Maxd  rowID
Maxd=size(out,2)-1;
%tic
flag=false; % if no intersection
dn=d;
outn=out;
Edgesn=Edges;




out_inter = lineSegmentIntersect(XY1,XY2);
Adj=out_inter.intAdjacencyMatrix;
[j1,j2]=find(Adj==1);


% further check if have shared nodes (that will give a,b=0,1) to avoid numerical error
%inter=0;
pair_inter=[];
for i=1:length(j1)
    a1=Edgesn(j1(i),1:2);
    a2=Edgesn(j2(i),1:2);
    if length(unique([a1,a2]))==4 && j1(i)<j2(i)
        pair_inter=[j1(i),j2(i),a1,a2];%[pair_inter;[j1(i),j2(i),a1,a2]];% edg j1(i) and edg j2(i) intersect edg1 [n1,n2] - edg2 [n3,n4] intersect
        break;
    end
end

%pair_inter
%% add new non-persistent node at intersection point
if ~isempty(pair_inter)
for i=1:1%size(pair_inter,1) %update 12_12_2016 only update one pair (as once updated the topology may change, intersection edges may change as well)
    %position
    ID=pair_inter(i,3:6);
    r=zeros(4,1);
    for j=1:4
        r(j)=rowID(ID(j));%find(out(:,1)==ID(j));
    end
    x=out_inter.intMatrixX(pair_inter(1,1),pair_inter(1,2));
    y=out_inter.intMatrixY(pair_inter(1,1),pair_inter(1,2));

    % calculate the distance of this [x,y] to four nodes in the edges
    % edg(j1(i),:),edg(j2(i),:) and define the closet one as the merging point
    D=pdist2([x,y],ps(r,1:2));
    D=[D',reshape(d(r),4,1)];
    [D,Ind]=sortrows(D,1);
    ix = find(D(:,2)<3, 1, 'first'); % merge to closest degree <=2 node
    if isempty(ix)
        ix=1;
    end
    nn=ID(Ind(ix));% merge at node nn

    % update links for node nn
    ou=zeros(1,Maxd+1);
    nei1=outn(r(Ind(ix)),2:end);
    nei=setdiff(union(nei1,ID),[nn,0]);
    ou(1:length(nei)+1)=[nn,nei];
    if length(ou)>Maxd+1
        Maxd=length(ou)-1;
        outn(:,end+1:Maxd+1)=0; %????
        %outn(:,1:Maxd+1)=0;
    end
    outn(r(Ind(ix)),:)=ou;
    dn(r(Ind(ix)))=length(nei);
    
    % update splitted edge nodes 
    if Ind(ix)<=2  %edg(j1(i),:), 
        % then split edge edg(j2(i),:) need to update r3,r4
        nei3=outn(r(3),2:end);
        nei=setdiff(union(nei3,nn),[0,ID(4)]); % ID(4)=out(r(4),1)
        ou=zeros(1,Maxd+1);
        ou(1:length(nei)+1)=[ID(3),nei];
        outn(r(3),:)=ou;
        dn(r(3))=length(nei);
        if d(r(3))==dn(r(3)) % no nn-n3 edge previously
            %XY1=[ps(r(Ind(ix)),1:2),ps(r(3),1:2)];
            %flag1=checkintersection(XY1,Lx,Ly);
            [len, flag1]=getDistanceCut(ps,r(Ind(ix)),r(3),Lx,Ly);
            Edgesn(end+1,:)=[sort([nn, ID(3)]),flag1,len];%norm(ps(r(Ind(ix)),1:2)-ps(r(3),1:2))];
        end
        
        nei4=outn(r(4),2:end);
        nei=setdiff(union(nei4,nn),[0,ID(3)]); % ID(4)=out(r(4),1)
        ou=zeros(1,Maxd+1);
        ou(1:length(nei)+1)=[ID(4),nei];
        outn(r(4),:)=ou;
        dn(r(4))=length(nei);
        if d(r(4))==dn(r(4)) % no nn-n4 edge previously
            %XY1=[ps(r(Ind(ix)),1:2),ps(r(4),1:2)];
            %flag1=checkintersection(XY1,Lx,Ly);
            [len, flag1]=getDistanceCut(ps,r(Ind(ix)),r(4),Lx,Ly);
            Edgesn(end+1,:)=[sort([nn, ID(4)]),flag1,len];%norm(ps(r(Ind(ix)),1:2)-ps(r(4),1:2))];
        end
        
        % update Edgesn after rewiring: change edge n3-n4 to n3-n1 and n4-n1 
        ind2=find(Edgesn(:,1)==ID(3) & Edgesn(:,2)==ID(4)); 
        Edgesn(ind2,:)=[];
        
        %distance
        
    elseif Ind(ix)>=3
        % then split edge edg(j1(i),:) need to update r1,r2
        nei1=outn(r(1),2:end);
        nei=setdiff(union(nei1,nn),[0,ID(2)]); % ID(4)=out(r(4),1)
        ou=zeros(1,Maxd+1);
        ou(1:length(nei)+1)=[ID(1),nei];
        outn(r(1),:)=ou;
        dn(r(1))=length(nei);
        if d(r(1))==dn(r(1)) % no nn-n1 edge previously
            %XY1=[ps(r(Ind(ix)),1:2),ps(r(1),1:2)];
            %flag1=checkintersection(XY1,Lx,Ly);
            [len, flag1]=getDistanceCut(ps,r(Ind(ix)),r(1),Lx,Ly);
            Edgesn(end+1,:)=[sort([nn, ID(1)]),flag1,len];%norm(ps(r(Ind(ix)),1:2)-ps(r(1),1:2))];
        end
        
        
        nei2=outn(r(2),2:end);
        nei=setdiff(union(nei2,nn),[0,ID(1)]); % ID(4)=out(r(4),1)
        ou=zeros(1,Maxd+1);
        ou(1:length(nei)+1)=[ID(2),nei];
        outn(r(2),:)=ou;
        dn(r(2))=length(nei); 
        if d(r(2))==dn(r(2)) % no nn-n2 edge previously
            %XY1=[ps(r(Ind(ix)),1:2),ps(r(2),1:2)];
            %flag1=checkintersection(XY1,Lx,Ly);
            [len, flag1]=getDistanceCut(ps,r(Ind(ix)),r(2),Lx,Ly);
            Edgesn(end+1,:)=[sort([nn, ID(2)]),flag1,len];%norm(ps(r(Ind(ix)),1:2)-ps(r(2),1:2))];
        end
        
        % update Edgesn after rewiring: change edge n3-n4 to n3-n1 and n4-n1
        ind2=find(Edgesn(:,1)==ID(1) & Edgesn(:,2)==ID(2)); 
        Edgesn(ind2,:)=[];
        
    end
  %  disp(strcat('intersection: merge edge ', num2str(ID(1)),'-',num2str(ID(2)),'and edge ',num2str(ID(3)),'-',num2str(ID(4)),'at node ',num2str(nn)));

end

%if ~isempty(pair_inter)
    flag=true;

end
%toc
