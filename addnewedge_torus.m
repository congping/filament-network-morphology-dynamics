function [pn,dn,outn,Ntn,Edgesn,flag]=addnewedge_torus(ps,d,out,Edges)
% add one new edge randomly
global Maxd  L rowID MaxID
Nt=size(ps,1);
Ntn=Nt;
pn=ps;
dn=d;
outn=out;
Edgesn=Edges;
flag=false;


Edges_temp=Edgesn;
Edges_temp(Edgesn(:,3)==1,:)=[];

S=Edges_temp(:,4);
s=cumsum(S);%s=cumsum(S(:,3));
r1=rand*s(end);
ix = find(s-r1>0,1);
n=Edges_temp(ix,[1,2,4]);%n=S(ix,:); % two neighbor nodes for the new edded node
n(4)=rowID(n(1));n(5)=rowID(n(2));
r2=rand;
t1=r2*ps(n(4),1:2)+(1-r2)*ps(n(5),1:2);
a=rand*2*pi;%random angle growing out from node t1
% initial set new edge of sufficent long
t2=t1+2.*L.*[cos(a),sin(a)];


% find the nearest intersection
XY1=[t1,t2];
XY2=[];%zeros(sum(d),4);
% list all edges and exclude duplicate edges
k=1;

for i=1:size(Edges_temp,1)
        ee1=rowID(Edges_temp(i,1));
        ee2=rowID(Edges_temp(i,2));
        XY2(k,:)=[ps(ee1,1:2),ps(ee2,1:2)];
        k=k+1;
end


out_inter = lineSegmentIntersect(XY1,XY2);
I= find(out_inter.intAdjacencyMatrix>0);
dd=out_inter.intNormalizedDistance1To2(I);
D=find(dd>0.0001);


if ~isempty(D) % update the position and neighbours and degree
    id=find(dd==min(dd(D)));
    ind=I(id(1));
    t2=[out_inter.intMatrixX(ind),out_inter.intMatrixY(ind)];
    %nn=edg(ind,:); % two neighbour node ID for the intersection point
    nn=Edges_temp(ind,1:2);
    % update the degree, neighbor, positions
    %new node position t1 at the edge n(1),n(2)
    Ntn=Ntn+1;% number of nodes
    %iNtn1=max(outn(:))+1;%new id
    iNtn1=MaxID+1;
    MaxID=MaxID+1;
    for m=1:dn(n(4))
        if outn(n(4),m+1)==n(2)
            outn(n(4),m+1)=iNtn1;
        end
    end
    for m=1:dn(n(5))
        if outn(n(5),m+1)==n(1)
            outn(n(5),m+1)=iNtn1;
        end
    end
    
    Ntn=Ntn+1;% number of nodes
    %iNtn2=max(outn(:))+1;%new id   
    iNtn2=MaxID+1;
    MaxID=MaxID+1;
    
    nn1=find(outn(:,1)==nn(1)); % give row number
    nn2=find(outn(:,1)==nn(2));
    for m=1:dn(nn1)
        if outn(nn1,m+1)==nn(2)
            outn(nn1,m+1)=iNtn2;
        end
    end
    for m=1:dn(nn2)
        if outn(nn2,m+1)==nn(1)
            outn(nn2,m+1)=iNtn2;
        end
    end
    
    
    dn(end+1:end+2)=[3,3];
    pn(end+1:end+2,1:3)=[[t1,iNtn1];[t2,iNtn2]];

    outn(end+1:end+2,:)=zeros(2,size(outn,2));
    outn(Ntn-1,1:4)=[iNtn1,n(1),n(2),iNtn2];
    outn(Ntn,1:4)=[iNtn2,nn(1),nn(2),iNtn1];
    %edg=[edg;[iNtn1,iNtn2]];
    [ind1,~,~]=checkCutBoundary(Edgesn,n(1),n(2));    
    [ind2,~,~]=checkCutBoundary(Edgesn,nn(1),nn(2));
    Edgesn(end+1,:)=[n(1),iNtn1,0,norm(t1-ps(n(4),1:2))];
    Edgesn(end+1,:)=[n(2),iNtn1,0,norm(t1-ps(n(5),1:2))];
    Edgesn(end+1,:)=[nn(1),iNtn2,0,norm(t2-ps(rowID(nn(1)),1:2))];
    Edgesn(end+1,:)=[nn(2),iNtn2,0,norm(t2-ps(rowID(nn(2)),1:2))];
    Edgesn(end+1,:)=[iNtn1,iNtn2,0,norm(t1-t2)];
    if length(ind1)~=1 || length(ind2)~=1
        disp('duplicated edges in addnewedg');
        ind1
        ind2
        pause
    else
        
    Edgesn([ind1,ind2],:)=[];
    flag=true;
    end
    display(strcat('add new edge ',num2str(iNtn1),' to ',num2str(iNtn2)));
    %checkGraph(pn,dn,outn,Edgesn);
end









