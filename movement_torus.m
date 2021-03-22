function [ps,Edges]=movement_torus(n,ps,d,out,Edges,Npp,Nt,dt,del,delta,stream,T_stream_on,T_stream_off,Lx,Ly)
global rowID
  for i=Npp+1:Nt %row
        %i
%         if i==2353
%             i
%         end
        dx=zeros(d(i),2);
        dx1=zeros(d(i),2);
        l=zeros(d(i),1);
        for j=2:d(i)+1
            j1=out(i,j);%ID
            jj1=find(out(:,1)==j1);
            %check if edge jj1-i cut boundary
            [~,~,ff]=checkCutBoundary(Edges,j1,out(i,1));
            newp=ps(jj1,1:2);
            if ff 
                if abs(ps(i,1)-newp(1))>Lx/2
                    if ps(i,1)>newp(1)
                        newp(1)=newp(1)+Lx;
                    else
                        newp(1)=newp(1)-Lx;
                    end
                end
                if abs(ps(i,2)-newp(2))>Ly/2
                    if ps(i,2)>newp(2)
                        newp(2)=newp(2)+Ly;
                    else
                        newp(2)=newp(2)-Ly;
                    end
                end                
            end
            dx(j-1,1:2)=newp-ps(i,1:2);%dx(j-1,:)=ps(jj1,1:2)-ps(i,1:2); 
            l(j-1)=norm(dx(j-1,1:2));
            
        end
        ind=find(l==0);
        dx(ind,:)=[];
        l(ind)=[];
        dx1=dx./l;
       % dx1(j-1,:)=dx(j-1,1:2)./l(j-1);%unit vector
        if length(l)==3% || d(i)==1 % only apply tension on the junction
             vel=sum(dx1,1); 
             psn=ps(i,1:2)+del*vel*dt;         
        elseif length(l)==1
             vel=sum(dx1,1); 
             psn=ps(i,1:2)+del*vel*dt;  
        elseif length(l)==2 % apply streaming in normal direction
                %     M= [((y3-y1)*y2^2+(y1^2-y3^2+(x1-x3)*(x1-2*x2+x3))*y2-y1^2*y3+(y3^2+(x2-x3)^2)*y1-y3*(x1-x2)^2)*((x1-x3)*y2+(x3-x2)*y1-y3*(x1-x2)),
                 %        -((y2-y3)*x1+(y3-y1)*x2+(y1-y2)*x3)*((x3-x1)*x2^2+(x1^2-x3^2+(y1-y3)*(y1-2*y2+y3))*x2-x3*x1^2+(x3^2+(y2-y3)^2)*x1-x3*(y1-y2)^2)];
                 %    M= [(-y1*y3^2+(x1^2+y1^2)*y3-x3^2*y1)*(x1*y3-x3*y1),
                 %        -(x1*y3-x3*y1)*(-x1*x3^2+(x1^2+y1^2)*x3-x1*y3^2)]% assuming P2=[0,0]
                 %M=-(-y1 . P3 . P3+P1 . P1 . y3,   x1 . P3 . P3-P1 . P1 . x3)
                     ccc=dx(1,1)*dx(2,2)-dx(2,1)*dx(1,2);
                     M=-ccc.*[-dx(1,2)*norm(dx(2,:))^2+norm(dx(1,:))^2*dx(2,2), dx(1,1)*norm(dx(2,:))^2-norm(dx(1,:))^2*dx(2,1)];
                      Kapp=2.*M./(norm(dx(1,:))^2.*norm(dx(2,:))^2.*norm(dx(1,:)-dx(2,:))^2);
                      % special case if points are on the same line
                  if abs(ccc)<1E-10
                                Kapp=[0,0];
                  end
                   if (n>=T_stream_on && n<=T_stream_off)
                    %if (ps(i,2)>v_stream_0 && ps(i,2)<v_stream_0+3) || (v_stream_0+3>L && ps(i,2)<mod(v_stream_0+3,L))        
                            %% calculate the streaming force for the valid component                           
                      if Kapp(1)==0 && Kapp(2)==0
                          sf=dx1(1,1)*stream(1)+dx1(1,2)*stream(2);%dot(dx1(1,:),stream);
                          N_stream=[dx1(1,2),-dx1(1,1)].*sign(sf);
                      else
                          N_stream=Kapp./norm(Kapp)*sign(Kapp(1));
                      end
                       temp=stream.*(ps(i,1)>0 && ps(i,1)<=Lx);%.*mean(l);% reset streming force to be proportional of the edge length link to the node.          
                     %ps(i,:)=ps(i,:)+v_stream.*dot(temp,v_stream).*dt;%stream.*dt;
                    else
                        temp=[0,0];%no streaming
                        N_stream=stream./norm(stream);
                    end
                    sf=temp(1)*N_stream(1)+temp(2)*N_stream(2);%dot(temp,N_stream);
                     psn=ps(i,1:2)+N_stream.*sf.*dt-delta.*Kapp.*dt;   
        elseif length(l)>3
             vel=sum(dx1); 
             psn=ps(i,1:2)+del*vel*dt;  
        end
        
        %check if any change of cutboundary due to node movement
       for dd=1:d(i)
           % E=sort([out(i,1),out(i,dd+1)]);
           % ind=find(Edges(:,1)==E(1) & Edges(:,2)==E(2));
            [ind,~,f1]=checkCutBoundary(Edges,out(i,1),out(i,dd+1));
            rr=rowID(out(i,dd+1));
            newp=ps(rr,1:2);
            if ~f1 % edge rowID i-rr within domain
               Edges(ind,4)=norm(psn-newp);
            else %if not within the domain update the postion for distance calculation                
                if abs(psn(1)-newp(1))>Lx/2
                if psn(1)>newp(1)
                    newp(1)=newp(1)+Lx;
                else
                    newp(1)=newp(1)-Lx;
                end
                end
                if abs(psn(2)-newp(2))>Ly/2
                if psn(2)>newp(2)
                    newp(2)=newp(2)+Ly;
                else
                    newp(2)=newp(2)-Ly;
                end    
                end
                Edges(ind,4)=norm(psn-newp);
            end
           
            if psn(1)<0 || psn(1)>Lx || psn(2)<0 || psn(2)>=Ly  %% move outsize the domain, may need to change the cut boundary, otherwise remain the same
                if ~f1
                    Edges(ind,3)=true;
                else % if edge i-rr cut bondary, now node i move out of the domain
                    newp(1)=mod(newp(1),Lx);
                    newp(2)=mod(newp(2),Ly);
                    psn_temp=psn;
                    if psn_temp(1)-newp(1)>Lx/2
                        psn_temp(1)=psn_temp(1)-Lx;
                    elseif newp(1)-psn_temp(1)>Lx/2
                        psn_temp(1)=psn_temp(1)+Lx;
                    end
                    if psn_temp(2)-newp(2)>Ly/2
                        psn_temp(2)=psn_temp(2)-Ly;
                    elseif newp(2)-psn_temp(2)>Ly/2
                        psn_temp(2)=psn_temp(2)+Ly;
                    end                    
                    XY1=[psn_temp,newp];%[ps([i,rr],1:2)];
                    Edges(ind,3)=checkintersection(XY1,Lx,Ly);
%                     XY2=[[0 0 Lx 0];[0 Ly 10 Ly]];
%                     outinterset = lineSegmentIntersect(XY1,XY2);
%                     dis=outinterset.intNormalizedDistance1To2;
%                     if (dis(1)>0 && dis(1)<1) || (dis(2)>0 && dis(2)<1)
%                         Edges(ind,3)=true;
%                     else
%                         Edges(ind,3)=false;
%                     end
                end

            end
       end   
        psn(1)=mod(psn(1),Lx); 
        psn(2)=mod(psn(2),Ly);
	    ps(i,1:2)=psn;    
  end
    
  
  %  checkGraph(ps,d,out,Edges);
        
        