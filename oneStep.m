function [ps,d,out,Nt,Edges]=oneStep(ps,d,out,Nt,Edges,Lx,Ly,Npp,Emin,del,delta,dt,n,T_crossoff, T_stream_on, T_stream_off, stream,rnew)
global rowID MaxID Maxd 

%% added nodes for curvature
[ps,d,out,Nt, Edges]=meshnodes_torus(ps,d,out,Nt, Edges,Lx,Ly,Emin);%[ps,d,out,Nt]=meshnodes_torus(ps,d,out,Nt, Emin);
rowID=getrowID(out);
    %    disp('afte meshnode')
        %cla(fig2);
        % plotnettorus(ps,out,Edges,Npp,Nt,Lx,false); 
%  v_stream_0=mod(v_stream_0+v_stream*dt,L);
   
    % now evolve Steiner points (Euler step) according to forces
% care of cut boundary check
   [ps,Edges]=movement_torus(n,ps,d,out,Edges,Npp,Nt,dt,del,delta,stream,T_stream_on,T_stream_off,Lx,Ly);
 %  disp('after movement')
       % cla(fig2);
        % plotnettorus(ps,out,Edges,Npp,Nt,Lx,false);     
    
    
    % then need to check to see if we need to split any new persistent points
    rr=rand;
    %flag0=false;
    if rr<rnew*dt  && n<T_crossoff % || n>=T_recross)%(n<T_stream_on|| n>=T_recross) %(n<T_stream_off || n>=T_recross)
        [pn,dn,outn,Ntn,Edgesn,flag0]=addnewedge_torus(ps,d,out,Edges);
        rowID=getrowID(outn);
        ps=pn;
        d=dn;
        out=outn;
        Nt=Ntn;
        Edges=Edgesn;
        [ps,d,out,Nt, Edges]=meshnodes_torus(ps,d,out,Nt, Edges,Lx,Ly,Emin);%[ps,d,out,Nt]=meshnodes_torus(ps,d,out,Nt, Emin);
        rowID=getrowID(out);
        %disp('try cross connection')
        %cla(fig2);
        % plotnettorus(ps,out,Edges,Npp,Nt,Lx,false); 
    end

       
       
    %% merge nodes
    [ps,d,out,Nt,Edges,flag1]=mergenodes_torus(ps,d,out,Nt,Lx,Ly,Edges);
    rowID=getrowID(out);
    %disp('after merging nodes')
   % cla(fig2);
   % plotnettorus(ps,out,Edges,Npp,Nt,Lx,false); 
       
       
    % check edge intersection
   [d,out,Edges,flag]=checkintersection_torus_fullregion(ps,d,out,Edges,Lx,Ly,Emin); % [d,out,flag]=checkintersection(ps,d,out);
    rowID=getrowID(out);
    flag2=flag;
    % after merge with intersection, may lead to new intersection
    mycount=0;
    while flag
       [d,out,Edges,flag]=checkintersection_torus_fullregion(ps,d,out,Edges,Lx,Ly,Emin); 
       rowID=getrowID(out);
       	mycount=mycount+1; 
    
        if mycount>30
            flag=false;
            disp('escaping...')
    %		size(ps)
        end
      %  disp('after checking intersection');
       % cla(fig2);
       % plotnettorus(ps,out,Edges,Npp,Nt,Lx,false); 
    end
    
    if  flag1 || flag2    %flag3 ||   
       % [ps,d,out,Edges]=split_torus(ps,d,out,Lx,Edges);
        % % identify nodes that need to have steiner points split off
        % (including degree >=4 non-persistent nodes and degree>=2 persistent nodes)
        [per,splits]=findsplits_torus(ps,d,out,Lx,Ly,Edges);
        %nsplit=1;
        while splits>0
            % split the nodes in "per" to create additional Steiner points
            [ps,d,out,Nt,Edges]=splitnodes_torus(ps,d,out,per,splits,Nt,Lx,Ly,Edges);
            rowID=getrowID(out);
            [per,splits]=findsplits_torus(ps,d,out,Lx,Ly,Edges);
           % nsplit=nsplit+1;
        %else
         %   ps=pp
        end
        %nsplit
      %  disp('after spliting');
      % cla(fig2);
      % plotnettorus(ps,out,Edges,Npp,Nt,Lx,false); 
    end