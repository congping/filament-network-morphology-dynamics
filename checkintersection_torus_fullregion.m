function [dn,outn,Edgesn,flag]=checkintersection_torus_fullregion(ps,d,out,Edges,Lx,Ly,Emin)
% this is checking potential edge intersection, if so, merge intersection
% point to neares node. So node number does not change, only the topolgy
% changes, out degree would increase for some nodes
global Maxd  rowID
rowID=getrowID(out);
Maxd=size(out,2)-1;
flag1=false;flag2=false;flag3=false;
%tic
%% check inner edges first
Edg=Edges;
loc=Edg(:,3)==0;
Edg(~loc,:)=[];

%k=1;
XY1=zeros(size(Edg,1),4);%[];%zeros(sum(d),4);
% list all edges and exclude duplicate edges
for i=1:size(Edg,1)
        ee1=rowID(Edg(i,1));
        ee2=rowID(Edg(i,2));
        XY1(i,:)=[ps(ee1,1:2),ps(ee2,1:2)];             
       % k=k+1;
end
XY2=XY1;

[dn,outn,Edg1,flag1]=checkintersection_torus(ps,d,out,Lx,Ly,Edg,XY1,XY2);
if flag1
Edgesn=[Edges(~loc,:);Edg1];
Edges=Edgesn;
rowID=getrowID(outn);
end



%% check left edges cut bounary
%% need to include all edges within [-3Emin, 3Emin]*[-3Emin, Ly+3Emin], y position of copied nodes depending on edges
pn=ps;
ind=find(pn(:,1)>Lx-2.5*Emin);
pn(ind,1)=pn(ind,1)-Lx;

Edg=[];
nodelist=pn(find(pn(:,1)<2.5*Emin),3); %this include nodes at [-3Emin,3Emin]
ind=false(size(Edges,1),1);
for i=1:size(Edges,1)
if ismember(Edges(i,1),nodelist) && ismember(Edges(i,2), nodelist) %isempty(setdiff(Edges(i,1:2),nodelist)) % only consider both end point of an edges is within the considered domain  
    ind(i)=true;%Edg=[Edg;Edges(i,:)];

end
end

Edg=Edges(ind,:);

if ~isempty(Edg)
XY1=[];%zeros(sum(d),4);
% list all edges and exclude duplicate edges
k=1;
for i=1:size(Edg,1)
        ee1=rowID(Edg(i,1));
        ee2=rowID(Edg(i,2));
        if pn(ee1,2)-pn(ee2,2) >Ly/2
             newp=pn(ee2,1:2);
             newp(2)=newp(2)+Ly;
             XY1(k,:)=[pn(ee1,1:2),newp];
        elseif pn(ee2,2)-pn(ee1,2)>Ly/2
             newp=pn(ee1,1:2);
             newp(2)=newp(2)+Ly;
             XY1(k,:)=[newp,pn(ee2,1:2)]; 
        else
            XY1(k,:)=[pn(ee1,1:2),pn(ee2,1:2)]; 
        end              
        k=k+1;
end
XY2=XY1;

[dn,outn,Edg2,flag2]=checkintersection_torus(pn,dn,outn,Lx,Ly,Edg,XY1,XY2);
if flag2
Edges=[Edges(~ind,:);Edg2];
rowID=getrowID(outn);
end
end


%% check bottom edges cut boundary
pn=ps;
ind=find(pn(:,2)>Ly-2.5*Emin);
pn(ind,2)=pn(ind,2)-Ly;

Edg=[];
nodelist=pn(find(pn(:,2)<2.5*Emin),3); %this include nodes at [-3Emin,3Emin]
ind=false(size(Edges,1),1);
for i=1:size(Edges,1)
if ismember(Edges(i,1),nodelist) && ismember(Edges(i,2), nodelist) %isempty(setdiff(Edges(i,1:2),nodelist)) % only consider both end point of an edges is within the considered domain  
    ind(i)=true;%Edg=[Edg;Edges(i,:)];

end
end

Edg=Edges(ind,:);
if ~isempty(Edg)
XY1=[];%zeros(sum(d),4);
% list all edges and exclude duplicate edges
k=1;
for i=1:size(Edg,1)
        ee1=rowID(Edg(i,1));
        ee2=rowID(Edg(i,2));
        if pn(ee1,1)-pn(ee2,1) >Lx/2
             newp=pn(ee1,1:2);
             newp(1)=newp(1)-Lx;
             XY1(k,:)=[newp,pn(ee2,1:2)]; 
        elseif pn(ee2,1)-pn(ee1,1)>Lx/2
             newp=pn(ee2,1:2);
             newp(1)=newp(1)-Lx;
             XY1(k,:)=[pn(ee1,1:2),newp];
        else
            XY1(k,:)=[pn(ee1,1:2),pn(ee2,1:2)]; 
        end              
        k=k+1;
end
XY2=XY1;


[dn,outn,Edg3,flag3]=checkintersection_torus(pn,dn,outn,Lx,Ly,Edg,XY1,XY2);
if flag3
Edges=[Edges(~ind,:);Edg3];
rowID=getrowID(outn);
end
  
end

Edgesn=Edges;  
flag=(flag1 || flag2 || flag3 );

%checkGraph(ps,dn,outn,Edgesn);
