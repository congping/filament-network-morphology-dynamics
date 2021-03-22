%% calculate ange of each ER filament
%%             90
%%             |
%%   180  -------------0
%%             |
%%            -90
                 
%%                                  H1
%% region R=[L1 L2,H1, H2]      L1   -----L2
%%                                  H2
function [angles,len]=getAngles_torus(p,Lx,Ly,Edges)
angles=[];
len=[];
% L1=R(1);
% L2=R(2);
% H1=R(3);
% H2=R(4);
global rowID
for i=1:size(Edges,1)
    p1=p(rowID(Edges(i,1)),1:2);
    p2=p(rowID(Edges(i,2)),1:2);
    if Edges(i,3)==1
        if p2(1)-p1(1)>Lx/2
            p1(1)=p1(1)+Lx;
        elseif p1(1)-p2(1)>Lx/2
            p1(1)=p1(1)-Lx;
        end
        if p2(2)-p1(2)>Ly/2
            p1(2)=p1(2)+Ly;
        elseif p1(2)-p2(2)>Ly/2
            p1(2)=p1(2)-Ly;
        end        
    end
    v=p1-p2;
    ang=atan2(v(2),v(1));% atan2(Y,X) in Four-quadrant inverse tangent
    if abs(ang)>pi/2
        ang=ang-sign(ang)*pi;%+2*pi;
    end
    angles=[angles,ang];
    len=[len,norm(v)];
end