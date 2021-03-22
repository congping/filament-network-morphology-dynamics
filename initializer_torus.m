function [pp,out,d, Edges]= initializer_torus(pp,Lx,Ly,flag)
pp(:,3)=1:size(pp,1);
pp0=pp;
Maxd=8;
Npp=size(pp,1);
bthresh=0.5;
zz=Npp+1;

% left and right
for qq=1:Npp    
    if pp(qq,1) < bthresh*Lx
        
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1)+Lx;
		pp(zz,2)=pp(qq,2);
		zz=zz+1;
    elseif  pp(qq,1) > (1-bthresh)*Lx
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1)-Lx;
		pp(zz,2)=pp(qq,2);
		zz=zz+1;
	end
end

% up and down
for qq=1:Npp    
    if pp(qq,2) < bthresh*Ly
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1);
		pp(zz,2)=pp(qq,2)+Ly;
		zz=zz+1;
    elseif pp(qq,2) > (1-bthresh)*Lx
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1);
		pp(zz,2)=pp(qq,2)-Ly;
		zz=zz+1;       
	end
end

% for corner
for qq=1:Npp    
    if pp(qq,2) < bthresh*Ly && pp(qq,1) < bthresh*Lx % region 1
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1)+Lx;
		pp(zz,2)=pp(qq,2)+Ly;
		zz=zz+1;
    elseif pp(qq,2) > (1-bthresh)*Ly && pp(qq,1) < bthresh*Lx  %region 3
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1)+Lx;
		pp(zz,2)=pp(qq,2)-Ly;
		zz=zz+1;     
    elseif pp(qq,2) < bthresh*Ly && pp(qq,1) > (1-bthresh) *Lx  %region 2
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1)-Lx;
		pp(zz,2)=pp(qq,2)+Ly;
		zz=zz+1;     
     elseif pp(qq,2) > (1-bthresh)*Ly && pp(qq,1) >(1-bthresh) *Lx  %region 4
        pp(zz,3)=zz;
		pp(zz,1)=pp(qq,1)-Lx;
		pp(zz,2)=pp(qq,2)+Ly;
		zz=zz+1; 
	end
end


Npp=size(pp,1);
d=zeros(Npp,1);
out=zeros(Npp,Maxd+1);
if flag %false %true
%pp=zeros(Npp,2);
w=zeros(Npp,Npp);
X=zeros(Npp,Npp);

% measure distances between them
for i=1:Npp
    for j=1:Npp
        if(i~=j)
        w(i,j)=sqrt((pp(i,1)-pp(j,1))^2+(pp(i,2)-pp(j,2))^2);
        X(i,j)=1;
        end
    end
end

%shortest distance between persistent nodes
%dmin=min(w(w>0));

% Kruskal's algorithm (code by Georgios Papachristoudis) finds
% minimal spanning tree: output are weights ws and adjacency matrix Xs
[~,~,Xs]=kruskal(X,w);

% collate node information for minimal spanning tree
out(1:Npp,1)=[1:Npp]';  %out = [ i, j1,j2,j3...]
%[i,j]=find(Xs>0)];
for i=1:Npp-1
    for j=i+1:Npp
        if(Xs(i,j)>0)
            % degree 
            d(i)=d(i)+1;
            d(j)=d(j)+1;
            % list of connections from node i/j
            out(i,d(i)+1)=j;
            out(j,d(j)+1)=i;
        end
    end
end



%% initialize with Delaunay triangulation
else
    TRI = delaunay(pp(:,1),pp(:,2));
    edges=[[TRI(:,1),TRI(:,2)];[TRI(:,1),TRI(:,3)];[TRI(:,3),TRI(:,2)]];
    edges=unique(edges,'rows');
    Adj=zeros(size(pp,1),size(pp,1));
    for i=1:size(edges,1)
        Adj(edges(i,1),edges(i,2))=1;
    end
    Adj=(Adj+Adj')>0;
    d=sum(Adj,2);
    for i=1:size(pp,1)
        j=find(Adj(i,:)>0);
        out(i,1:length(j)+1)=[i,j];
    end
    
    
end




%% back to the network in the domain only

Edges=outToEdges(out);
Edges=sort(Edges,2);
Edges=unique(sort(Edges,2),'rows');

%% identify cut of boundary: cut along x, cut along y, cut along both x and y
Npp=size(pp0,1);
ind1=find(Edges(:,1)<=Npp & Edges(:,2)>Npp);  % if one node out of domain
Edges(ind1,3)=true;
Edges(ind1,2)=mod((Edges(ind1,2)-1),Npp)+1;

ind2=find(Edges(:,1)>Npp & Edges(:,2)>Npp);
%check position
for i=1:length(ind2)
    nn=Edges(ind2(i),1:2);
    f1=(pp(nn,1)<0) .* ((pp(nn,2)>=0) .* (pp(nn,2)<Ly));% in left outbord
    f2=(pp(nn,1)>=Lx) .* ((pp(nn,2)>=0) .* (pp(nn,2)<Ly));% in right outbord
    f3=(pp(nn,2)>=Ly) .* ((pp(nn,1)>=0) .* (pp(nn,1)<Lx));% in up outbord
    f4=(pp(nn,2)<0)  .* ((pp(nn,1)>=0) .* (pp(nn,1)<Lx));% in  bottom outbord
    f=[f1,f2,f3,f4];
    if sum(f(:))==2 && sum(sum(f)<2)==4
        XY1=[pp(nn(1),1:2),pp(nn(2),1:2)];
        ff=checkintersection(XY1,Lx,Ly);
        Edges(ind2(i),3)=ff;
        if ff
             Edges(ind2(i),1:2)=mod((Edges(ind2(i),1:2)-1),Npp)+1;
        end
%         XY2=[[0 0 Lx 0];[0 Ly Lx Ly]];
%         outinterset = lineSegmentIntersect(XY1,XY2);
%         dis=outinterset.intNormalizedDistance2To1;
%         if sum((dis<1).* (dis>0))
%                 Edges(ind2(i),3)=true;
%                 Edges(ind2(i),1:2)=mod((Edges(ind2(i),1:2)-1),Npp)+1;
%         end
    end    
end
ind3=find(Edges(:,1)>Npp & Edges(:,2)>Npp);
Edges(ind3,:)=[];



b=sort(Edges(:,1:2),2);
[~,ic]=unique(b,'rows');
Edges=Edges(ic,:);
Edges(:,1:2)=sort(Edges(:,1:2),2);
[out,d]=EdgesToOut(Npp,Edges);

global rowID
rowID=getrowID(out);
for i=1:size(Edges,1)
    nn=Edges(i,1:2);%node id
    len=getDistance(pp,Edges,nn,Lx,Ly);
    Edges(i,4)=len;
end


pp=pp0;


