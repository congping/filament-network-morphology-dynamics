function [pn,dn,outn,Ntn,Edgesn]=splitnodes_torus(ps,d,out,per,splits,Nt,Lx,Ly,Edges)
% split the nodes marked in per to create new network
% ONLY split persistent nodes
% ps- persistent node position
%per triple element for split [i,j,k] angle [i,j]-[i,k] 
global Maxd eps rowID MaxID
Maxd=size(out,2)-1;
Ntn=Nt+splits;
%pn=zeros(Ntn,3);
%dn=zeros(Ntn,1);
%outn=zeros(Ntn,Maxd+1);
Edgesn=Edges;
pn(1:Nt,:)=ps(1:Nt,:);
dn(1:Nt)=d(1:Nt);
outn(1:Nt,:)=out(1:Nt,:);

N=Nt;
% for each point to be split, add a new steiner point
for k=1:splits
    iN=MaxID+1;
    MaxID=MaxID+1;
    ii=per(k,1); %%ID
    jj1=per(k,2);%%ID
    jj2=per(k,3);%%ID
    dx1=per(k,7:8);
    dx2=per(k,9:10);
    i=find(outn(:,1)==ii);%% row number
    j1=find(outn(:,1)==jj1);%% row number
    j2=find(outn(:,1)==jj2);%% row number
    %i=rowID(ii);j1=rowID(jj1);j2=rowID(jj2);
    % create an extra node
    N=N+1;
    
    % create new point close to ps(i,:)
    pn(N,1:2)=pn(i,1:2)+eps.*(dx1+dx2)./norm(dx1+dx2);
    pn(N,3)=iN;
    %new edges iN-ii, iN-jj1, iN-jj2, delete edge ii-jj1, ii-jj2
    [indj1,~,cutf1]=checkCutBoundary(Edgesn,ii,jj1);
    [indj2,~,cutf2]=checkCutBoundary(Edgesn,ii,jj2);
    newp1=pn(j1,1:2);
    newp2=pn(j2,1:2);
    if cutf1
        if abs(pn(i,1)-pn(j1,1))>Lx/2
            if pn(i,1)>pn(j1,1)  
                newp1(1)=newp1(1)+Lx;
            else
                newp1(1)=newp1(1)-Lx;
            end
        end
        if abs(pn(i,2)-pn(j1,2))>Lx/2
            if pn(i,2)>pn(j1,2)  
                newp1(2)=newp1(2)+Ly;
            else
                newp1(2)=newp1(2)-Ly;
            end
        end
    end
    dj1=norm(pn(N,1:2)-newp1);
        
    if cutf2
        if abs(pn(i,1)-pn(j2,1))>Lx/2
            if pn(i,1)>pn(j2,1)  
                newp2(1)=newp2(1)+Lx;
            else
                newp2(1)=newp2(1)-Lx;
            end
        end
        if abs(pn(i,2)-pn(j2,2))>Lx/2
            if pn(i,2)>pn(j2,2)  
                newp2(2)=newp2(2)+Ly;
            else
                newp2(2)=newp2(2)-Ly;
            end
        end
    end
    dj2=norm(pn(N,1:2)-newp2);
    di=norm(pn(N,1:2)-pn(i,1:2));
    
   
    if pn(N,1)>=0 && pn(N,1)<Lx && pn(N,2)>=0 && pn(N,2)<Ly %% within the domain
        Edgesn(end+1,:)=[ii,iN,0,di];
        Edgesn(end+1,:)=[jj1,iN,cutf1,dj1];
        Edgesn(end+1,:)=[jj2,iN,cutf2,dj2];    
    else
        Edgesn(end+1,:)=[ii, iN, 1, di];
        
        if ~cutf1
            Edgesn(end+1,:)=[jj1,iN,~cutf1,dj1];
        else
            XY1=[pn(i,1:2),pn(j1,1:2)];%[pn([i,j1],1:2)];
            ff=checkintersection(XY1,Lx,Ly);
%             XY2=[[0 0 Lx 0];[0 Ly 10 Ly]];
%             outinterset = lineSegmentIntersect(XY1,XY2);
%             dis=outinterset.intNormalizedDistance2To1;
            if ff
                Edgesn(end+1,:)=[jj1,iN,true,dj1];
            else
                Edgesn(end+1,:)=[jj1,iN,false,dj1];
            end
        end     
        
             
        if ~cutf2
            Edgesn(end+1,:)=[jj2,iN,~cutf2,dj2];
        else
            XY1=[pn(i,1:2),pn(j2,1:2)];%[pn([i,j2],1:2)];
            ff=checkintersection(XY1,Lx,Ly);
%             XY2=[[0 0 Lx 0];[0 Ly 10 Ly]];
%             outinterset = lineSegmentIntersect(XY1,XY2);
%             dis=outinterset.intNormalizedDistance2To1;
            if ff%dis>0 && dis<1
                Edgesn(end+1,:)=[jj2,iN,true,dj2];
            else
                Edgesn(end+1,:)=[jj2,iN,false,dj2];
            end
        end   
        
        pn(N,1)=mod(pn(N,1),Lx);
        pn(N,2)=mod(pn(N,2),Ly);
    end
    Edgesn([indj1,indj2],:)=[];
    
    
    % repack to remove connections i to j1 and j2
    j=1;  
    % changed k to kk by Lin 10-10-2016
    for kk=1:dn(i)
        if ((outn(i,kk+1)~=jj1) && (outn(i,kk+1)~=jj2))
            outn(i,j+1)=outn(i,kk+1);
            j=j+1;
        end
    end
    % one less connection
    outn(i,dn(i)+1)=0;
    dn(i)=dn(i)-1;
    % add connection i to New ID  iN
    %iN=max(outn(:,1))+1;
    
    outn(i,dn(i)+1)=iN;

    % replace connections from j1,j2 -> i with j1,j2 -> N
     
    % changed k-kk by Lin 10-10-2016
    for kk=1:dn(j1)
        if outn(j1,kk+1)==ii
            outn(j1,kk+1)=iN;
        end
    end
    for kk=1:dn(j2)
        if outn(j2,kk+1)==ii
            outn(j2,kk+1)=iN;
        end
    end
    
    % new connections from N to i, j1 and j2
    dn(N)=3;
    outn(N,1)=iN;% new ID
    outn(N,1+1)=ii;
    outn(N,2+1)=jj1;
    outn(N,3+1)=jj2;
    
%     % renumber any remaining common nodes to those that must also be split
%updated accordingly
    for l=k+1:splits
        for m=2:3
            if (per(l,1)==jj1 && per(l,m)==ii) || (per(l,1)==jj2 && per(l,m)==ii)
                per(l,m)=iN;per(l,m+3)=N;
            end
        end
    end
    
end
end