function [ps,d,out,Edges,Nt,Npp]=inputgraph(graph,Npp,k)
ps=graph(k).position;
d=graph(k).degree;
out=graph(k).outnode;
Edges=graph(k).edges;
Nt=size(ps,1);
