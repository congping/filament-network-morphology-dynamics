function [len, flag]=getDistanceCut(ps,r1,r2,Lx,Ly)
 psnew=ps(r2,1:2);
% cut cross x
if ps(r2,1)-ps(r1,1)>0.5*Lx
    psnew(1)=psnew(1)-Lx;
elseif ps(r1,1)-ps(r2,1)>0.5*Lx        
    psnew(1)=psnew(1)+Lx;
end  

% cut cross y
if ps(r2,2)-ps(r1,2)>0.5*Ly
    psnew(2)=psnew(2)-Ly;
elseif ps(r1,2)-ps(r2,2)>0.5*Ly        
    psnew(2)=psnew(2)+Ly;
end    
len = norm(ps(r1,1:2)-psnew);  

XY1=[ps(r1,1:2),psnew];
flag=checkintersection(XY1,Lx,Ly);