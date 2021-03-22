function plotnet(p,out,Edges,Np,flagnode)
cla;
hold on;
global rowID
rowID=getrowID(out);
for ii=1:size(p,1)
        if flagnode
       %text(p(ii,1),p(ii,2),num2str(rowID(out(ii,1))),'FontSize',10);
        text(p(ii,1),p(ii,2),num2str(p(ii,3)),'FontSize',10);
        hold on
        end
end
for i=1:size(Edges,1)
    %getrowId
    ii=Edges(i,1);%node ID %% need to return to row ID
    jj=Edges(i,2);%node ID 
    ri=rowID(ii);
    rj=rowID(jj);
    plot([p(ri,1),p(rj,1)],[p(ri,2),p(rj,2)],'b','LineWidth',2);
         %hold on
   
end


plot(p(1:Np,1),p(1:Np,2),'r.','MarkerSize',8);

hold off;
box on
drawnow;

