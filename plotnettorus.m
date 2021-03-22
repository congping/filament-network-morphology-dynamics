function plotnettorus(p,out,Edges,Np,Lx,Ly,flagnode,color)
%cla;
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
    if Edges(i,3)==0
         plot([p(ri,1),p(rj,1)],[p(ri,2),p(rj,2)],'color',color,'LineWidth',2);
         %hold on
    else       
        %% cross x        
        if abs(p(ri,1)-p(rj,1))>Lx*0.5 && abs(p(rj,2)-p(ri,2))<Ly*0.5
            if p(rj,1)>p(ri,1)
                    plot([p(ri,1)+Lx,p(rj,1)],[p(ri,2),p(rj,2)],'color',color,'LineWidth',2);
                    plot([p(ri,1),p(rj,1)-Lx],[p(ri,2),p(rj,2)],'color',color,'LineWidth',2);
            elseif p(rj,1)<p(ri,1)
                    plot([p(ri,1)-Lx,p(rj,1)],[p(ri,2),p(rj,2)],'color',color,'LineWidth',2);
                    plot([p(ri,1),p(rj,1)+Lx],[p(ri,2),p(rj,2)],'color',color,'LineWidth',2);
            end
            
            %% cross y
        elseif abs(p(ri,1)-p(rj,1))<Lx*0.5 && abs(p(rj,2)-p(ri,2))>Ly*0.5
            if p(rj,2)>p(ri,2)
                    plot([p(ri,1),p(rj,1)],[p(ri,2),p(rj,2)-Ly],'color',color,'LineWidth',2);
                    plot([p(ri,1),p(rj,1)],[p(ri,2)+Ly,p(rj,2)],'color',color,'LineWidth',2);
            elseif p(rj,2)<p(ri,2)
                    plot([p(ri,1),p(rj,1)],[p(ri,2),p(rj,2)+Ly],'color',color,'LineWidth',2);
                    plot([p(ri,1),p(rj,1)],[p(ri,2)-Ly,p(rj,2)],'color',color,'LineWidth',2);
            end  
          
              %% cross both
        elseif abs(p(ri,1)-p(rj,1))>Lx*0.5 && abs(p(rj,2)-p(ri,2))>Ly*0.5
            if (p(rj,2)<p(ri,2) && p(rj,1)<p(ri,1)) || (p(rj,2)>p(ri,2) && p(rj,1)>p(ri,1))
                    p1=[min(p(ri,1),p(rj,1)),min(p(ri,2),p(rj,2))];
                    p2=[max(p(ri,1),p(rj,1)),max(p(ri,2),p(rj,2))];                    
                    plot([p1(1)+Lx,p2(1)],[p1(2)+Ly,p2(2)],'color',color,'LineWidth',2);
                    plot([p1(1)+Lx,p2(1)],[p1(2),p2(2)-Ly],'color',color,'LineWidth',2);
                    plot([p1(1),p2(1)-Lx],[p1(2)+Ly,p2(2)],'color',color,'LineWidth',2);
                    plot([p1(1),p2(1)-Lx],[p1(2),p2(2)-Ly],'color',color,'LineWidth',2);
            elseif (p(rj,2)>p(ri,2) && p(rj,1)<p(ri,1)) || (p(rj,2)<p(ri,2) && p(rj,1)>p(ri,1))
                    p1=[min(p(ri,1),p(rj,1)),max(p(ri,2),p(rj,2))];
                    p2=[max(p(ri,1),p(rj,1)),min(p(ri,2),p(rj,2))];
                    plot([p1(1),p2(1)-Lx],[p1(2),p2(2)+Ly],'color',color,'LineWidth',2);
                    plot([p1(1)+Lx,p2(1)],[p1(2)-Ly,p2(2)],'color',color,'LineWidth',2);
                    plot([p1(1),p2(1)-Lx],[p1(2)-Ly,p2(2)],'color',color,'LineWidth',2);
                    plot([p1(1)+Lx,p2(1)],[p1(2),p2(2)+Ly],'color',color,'LineWidth',2);
            else
                disp('plotting error!')
                pause;
            end             
        end

    end
end


plot(p(1:Np,1),p(1:Np,2),'r.','MarkerSize',8);

hold off;
box on
drawnow;

