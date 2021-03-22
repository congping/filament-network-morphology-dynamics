%% edgeflag true-is an edge
%% cutflag ture - cut boundary
%% ind  - the row index of the edge in the list Edges


function [ind,edgeflag,cutflag]=checkCutBoundary(Edges,e1,e2)
E=sort([e1,e2]);

ind=find(Edges(:,1)==E(1) & Edges(:,2)==E(2));
%loc= Edges(:,1)==E(1);
%ind1=find(Edges(loc,2)==E(2));
%ind=loc(ind1);
if length(ind)~=1
    edgeflag=false;
    cutflag=false;
else
    edgeflag=true;
    cutflag=Edges(ind,3);
end
    
