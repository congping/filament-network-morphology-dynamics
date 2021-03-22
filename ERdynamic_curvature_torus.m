%% this matlab code run network morphology dynamics underlying stream

clear all
tic
global Maxd eps eta L Npp rowID MaxID 
flagplot=true; %false;
% domain of the torus
Lx=20;%
Ly=20;%
L=max(Lx,Ly);

% set time step and run time and stream on and off time and recorss on and off time (if there is any)
Tmax=10;
dt=0.0002;%0.02;
%T_recross=100;%260;
T_stream_off=4;%200;%160;
T_stream_on=1;%0;%150;%80;
T_crossoff=0;%0;


% no of persistent points (anchors)
Npp=150;
MaxID=Npp; % MaxID is for node, maximal ID has been used in the evolution

% max degree
Maxd=8;

%calculated variables
prop=[];%proportion of weighted angle paralle to stream

% generate uniform distribution of persistent points (anchors)
 for i=1:Npp
    pp(i,1)=rand*Lx;
    pp(i,2)=rand*Ly;
 end

%% if true initialize the graph using minimal spanning tree
[ps,out,d,Edges]= initializer_torus(pp,Lx,Ly,true);
Nt=size(ps,1); % number of nodes which varies in the evolution

%dmin=pdist(pp);
dmin=min(Edges(:,4));
MaxID=max(out(:,1));

%eta distance to merge 
eta=0.2*min(0.1,dmin);%0.1*min(0.12,dmin);%0.005*L;

% eps to shift (in splitting)
eps=2*eta;%0.008;

% del tension and drag
del=1;%b=1; %F/ (6pi eta R)

% del*dt  < eta < eps
delta=1;% tan/sigma
% new edge creation rate
rnew = 0;%2*0.0048*Lx*Ly;

stream=5.*[1,0];%[cos(pi/3),sin(pi/3)];%[1,0]; streaming velocity


Emin = 0.4; % for adding extra nodes for curvature  Emin<pi tention / streaming velocity / drag coeff

while max(Edges(:,4))>Emin
[ps,d,out,Nt, Edges]=meshnodes_torus(ps,d,out,Nt, Edges,Lx,Ly,Emin);
end

%% plot orginal graph
fig2=figure;
clf;
plotnettorus(ps,out,Edges,Npp,Lx,Ly,false,[0 1 0])
%title('Initial graph')
axis equal
axis([0 Lx 0 Ly])
set(fig2,'Parent');


%% initial split
rowID=getrowID(out);
[ps,d,out,Edges]=split_torus(ps,d,out,Lx,Ly,Edges);
Nt=size(out,1);

% main loop for timestepping
kk=1;
for n=0:dt:Tmax
      n
    
    [ps,d,out,Nt,Edges]=oneStep(ps,d,out,Nt,Edges,Lx,Ly,Npp,Emin,del,delta,dt,n,T_crossoff, T_stream_on, T_stream_off, stream,rnew);
     if n==T_stream_off
         n
     off_graph.position=ps;
     off_graph.degree=d;
     off_graph.outnode=out;
     off_graph.edges=Edges;
     off_graph.Nt=Nt;
     end

    
    if mod(kk,500)==1 
        [angles,len]=getAngles_torus(ps,Lx,Ly,Edges);
        ind=abs(angles-atan2(stream(2),stream(1)))<pi/4;
        prop=[prop;[sum(len(ind)),sum(len(~ind)),sum(len),sum(len(ind))/sum(len),n]];
    end
    
        % display current configuration
    if flagplot  && mod(kk,500)==1 %&& n>100
     %  figure(2);
       cla(fig2);
       plotnettorus(ps,out,Edges,Npp,Lx,Ly,false,[0 1 0]); 
       s=strcat('t=',num2str(n),' s');
       text(0.15,1,s,'FontSize',14);
        F(floor(kk/500)+1)=getframe;
    end
    kk=kk+1;

end


    endgraph.position=ps;
    endgraph.degree=d;
    endgraph.outnode=out;
    endgraph.edges=Edges;
    endgraph.Nt=Nt;
str=strcat('angles_len_L_',num2str(L),'Npp',num2str(Npp),'_dt',num2str(dt),'_stream',num2str(stream(1)),'_',num2str(stream(2)),'_eta',num2str(eta),'_del',num2str(del),'_delta',num2str(delta),'_eps',num2str(eps),'_Emin_',num2str(Emin),'r_cross',num2str(rnew),'_T_crossoff',num2str(T_crossoff),'T_stream_on',num2str(T_stream_on),'_T_stream_off',num2str(T_stream_off),'_Tmax',num2str(Tmax),'.mat');

save(str,'prop','ini_graph','endgraph','off_graph');
if flagplot
    str=strcat('movie',str(11:end));
    v = VideoWriter(str);
    open(v)
    writeVideo(v,F);
    close(v)

end

toc