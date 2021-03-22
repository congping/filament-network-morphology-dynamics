
function [per,splits]=findsplits_torus(ps,d,out,Lx,Ly,Edges)
% Find nodes that need to be split because of an angle < 2pi/3
% returns either per(1,1)=0 no split, or need to split
% edges (per(i,1),per(i,2)) and (per(i,1),per(i,3)).
% splits= number of nodes need to split
% pp -  position of persistent nodes
global Npp
%tic
ind1=find(d(1:Npp)>=2);% find split for persistent nodes (degree>=2)

ind2=Npp+find(d(Npp+1:end)>=4);% find split for non-persistent nodes (require degree >=4)

ind=union(ind1,ind2);

splits=0;
per=[];
for ii=1:length(ind)
    %ii
    ri=ind(ii);
    [split,oneper]=find_onesplit_torus(ps,d,out,ri,Lx,Ly,Edges);
    if split>0
       splits=splits+split;
       per=[per;oneper];
    end
end