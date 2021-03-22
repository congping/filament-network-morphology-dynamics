function [pn,dn,outn,Ntn, Edgesn]=meshnodes_torus(ps,d,out,Nt, Edges,Lx,Ly,Emin)

global MaxID 

% division threshold
thresh=Emin; 

% start with old
pn=ps;
dn=d;
outn=out;
Ntn=Nt;
Edgesn=Edges;
%maxd=8;
%j=4;
global rowID
rowID=getrowID(out);
MaxID=max(out(:,1));

ind=find(Edges(:,4)>thresh);


for i=1:length(ind)
    MaxID=MaxID+1;
    Ntn=Ntn+1;
    id1=Edges(ind(i),1);
    id2=Edges(ind(i),2);
    r1=rowID(id1);
    r2=rowID(id2);
    k=find(out(r1,:)==id2);
    outn(r1,k)=MaxID; % connection from to new
    ii=find(outn(r2,:)==id1); % which connection was coming back?
    outn(r2,ii)=MaxID;
    %out(r2,ii)=0; % don't do this again...
    dn(Ntn)=2;
    outn(Ntn,1)=MaxID;
    outn(Ntn,2)=id1;
    outn(Ntn,3)=id2;
    pn(Ntn,3)=MaxID;
    if Edges(ind(i),3)==0
            pn(Ntn,1:2)=(ps(r1,1:2)+ps(r2,1:2))./2;
            Edgesn(end+1,:) = [id1, MaxID, 0, norm(pn(Ntn,1:2)-ps(r1,1:2))];
            Edgesn(end+1,:) = [id2, MaxID, 0, norm(pn(Ntn,1:2)-ps(r2,1:2))];
    elseif Edges(ind(i),3)==1      
           %rewiring ps(r2,:) to the torus
           newps=ps(r2,1:2);
           if abs(ps(r1,1)-ps(r2,1)) > Lx*0.5
           if ps(r1,1)>ps(r2,1)
               newps(1)=newps(1)+Lx;
           else
               newps(1)=newps(1)-Lx;
           end
           end
           if abs(ps(r1,2)-ps(r2,2))>Ly*0.5
           if ps(r1,2)>ps(r2,2)
               newps(2)=newps(2)+Ly;
           else
               newps(2)=newps(2)-Ly;
           end
           end
           pn(Ntn,1:2)=(ps(r1,1:2)+newps)./2;
           
           if (pn(Ntn,1)>=0 && pn(Ntn,1)<Lx)   && ( pn(Ntn,2)>=0 && pn(Ntn,2)<Ly)    %% new point is within the domain     
                Edgesn(end+1,:) = [id1, MaxID, 0, norm(pn(Ntn,1:2)-ps(r1,1:2))];
                Edgesn(end+1,:) = [id2, MaxID, 1, norm(pn(Ntn,1:2)-newps)];
           else           
                Edgesn(end+1,:) = [id1, MaxID, 1, norm(pn(Ntn,1:2)-ps(r1,1:2))];
                %ff=(pn(Ntn,1)>=0 && pn(Ntn,1)<Lx) || (pn(Ntn,2)>=0 && pn(Ntn,2)<Ly);
                dd=norm(pn(Ntn,1:2)-newps);

               % fixed r2
                if pn(Ntn,1)-ps(r2,1)>Lx/2
                    pn(Ntn,1)=pn(Ntn,1)-Lx;
                elseif pn(r2,1)-pn(Ntn,1)>Lx/2
                    pn(Ntn,1)=pn(Ntn,1)+Lx;
                end
                if pn(Ntn,2)-ps(r2,2)>Ly/2
                    pn(Ntn,2)=pn(Ntn,2)-Ly;
                elseif pn(r2,2)-pn(Ntn,2)>Ly/2
                    pn(Ntn,2)=pn(Ntn,2)+Ly;
                end            
                 XY1=[pn(Ntn,1:2),ps(r2,1:2)];
                 ff=checkintersection(XY1,Lx,Ly);
                %XY2=[[0 0 Lx 0];[0 Ly 10 Ly]];
                %outinterset = lineSegmentIntersect(XY1,XY2);
                %dis=outinterset.intNormalizedDistance1To2;
                %if (dis(1)>0 && dis(1)<1) || (dis(2)>0 && dis(2)<1)
                 %   ff=true;
                %else
                 %   ff=false;
                %end
                 Edgesn(end+1,:) = [id2, MaxID, ff,dd ];
                 pn(Ntn,1)=mod(pn(Ntn,1),Lx);
                 pn(Ntn,2)=mod(pn(Ntn,2),Ly);
           end            
    end
end
Edgesn(ind,:)=[];


%checkGraph(pn,dn,outn,Edgesn);