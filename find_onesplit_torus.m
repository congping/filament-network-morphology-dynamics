function  [splits,per]=find_onesplit_torus(ps,d,out,ri,Lx,Ly,Edges)
global rowID
    % find split for one particular node ID IDi row ri
IDi=out(ri,1);
splits=0;    
per=zeros(1,6);    
jj=zeros(d(ri),1);
dp=zeros(d(ri),2);
ndp=zeros(d(ri),2);
for kk=1:d(ri)
    jj(kk)=out(ri,kk+1);  %ID
    ix=rowID(jj(kk));
    %check if edge IDi jj cuts boundary
    [~,edgeflag,cutflag]=checkCutBoundary(Edges,IDi,jj(kk));
    if ~edgeflag
        disp('not an edge in find_onesplit_periodic!');
        break
    elseif cutflag
        %cut boundary
        newpx=ps(ix,1:2);
        if abs(ps(ri,1)-ps(ix,1))>Lx*0.5
            if ps(ri,1)>ps(ix,1)             
                newpx(1)=newpx(1)+Lx;
            else
                newpx(1)=newpx(1)-Lx;
            end
        elseif abs(ps(ri,2)-ps(ix,2))>Ly*0.5
            if ps(ri,2)>ps(ix,2)             
                newpx(2)=newpx(2)+Ly;
            else
                newpx(2)=newpx(2)-Ly;
            end  
        end
        dp(kk,:)=newpx-ps(ri,1:2);
    elseif ~cutflag
        dp(kk,:)=ps(ix,1:2)-ps(ri,1:2);
    end
    
    ndp(kk,:)=dp(kk,1:2)/norm(dp(kk,1:2));
end

maxca=-1;

for k=1:d(ri)-1
    for l=k+1:d(ri)
        % cosine of separating angle
        ca =ndp(k,:)*ndp(l,:)';
        if ca>maxca
            maxca=ca;
            ikk=jj(k);%ID
            ill=jj(l);%ID
            rkk=rowID(jj(k));
            rll=rowID(jj(l));
            dp1=dp(k,:);
            dp2=dp(l,:);
            %rkk=kk;
            %rll=ll;
        end
    end
end

if maxca>cos(2*pi/3)+0.1
        % if maxca is big enough then mark edges
        % (i,mkk) and (i,mll)) to be split.
        splits=splits+1;
        per=[IDi,ikk,ill,ri,rkk,rll, dp1,dp2]; % IDi,IDj1 IDj2, ri,rj1,rj2, vector1,vector2]
end
        