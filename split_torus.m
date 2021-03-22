%% initial split
% identify persistent nodes that need to have steiner points split off
function [ps,d,out,Edges]=split_torus(ps,d,out,Lx,Ly,Edges)
global rowID
[per,splits]=findsplits_torus(ps,d,out,Lx,Ly,Edges);
Nt=size(ps,1);

while splits>0    
    % split the nodes in "per" to create additional Steiner points
    [ps,d,out,Nt,Edges]=splitnodes_torus(ps,d,out,per,splits,Nt,Lx,Ly,Edges);
    rowID=getrowID(out);
    [per,splits]=findsplits_torus(ps,d,out,Lx,Ly,Edges);
   % checkGraph(ps,d,out,Edges);
end


