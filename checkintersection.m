function flag=checkintersection(XY1,Lx,Ly)
% XY1=[pn(r1,1:2),newp];%[ps([i,rr],1:2)];
XY2=[[0 0 Lx 0];[0 Ly 10 Ly]];
outinterset = lineSegmentIntersect_1(XY1,XY2);
dis=outinterset.intNormalizedDistance1To2;
if (dis(1)>0 && dis(1)<1) || (dis(2)>0 && dis(2)<1)
    flag1=true;
else
    flag1=false;
end

XY2=[[0 0 0 Ly];[Lx 0 Lx Ly]];
outinterset = lineSegmentIntersect_1(XY1,XY2);
dis=outinterset.intNormalizedDistance1To2;
if (dis(1)>0 && dis(1)<1) || (dis(2)>0 && dis(2)<1)
    flag2=true;
else
    flag2=false;
end

flag=flag1 || flag2;