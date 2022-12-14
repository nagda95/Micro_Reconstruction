function Vol=Volume_Geode(R,Structure,Angles)
Vol=0;
for i=1:size(Structure,1)
    Rmoy=(R(Structure(i,1),1)+R(Structure(i,2),1)+R(Structure(i,3),1))/3;
    phi1=Angles(Structure(i,1),2);
    teta1=Angles(Structure(i,1),1);
    x1=Rmoy*sin(phi1)*cos(teta1);
    y1=Rmoy*sin(phi1)*sin(teta1);
    z1=Rmoy*cos(phi1);
    phi2=Angles(Structure(i,2),2);
    teta2=Angles(Structure(i,2),1);
    x2=Rmoy*sin(phi2)*cos(teta2);
    y2=Rmoy*sin(phi2)*sin(teta2);
    z2=Rmoy*cos(phi2);
    phi3=Angles(Structure(i,3),2);
    teta3=Angles(Structure(i,3),1);
    x3=Rmoy*sin(phi3)*cos(teta3);
    y3=Rmoy*sin(phi3)*sin(teta3);
    z3=Rmoy*cos(phi3);
    a=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    b=sqrt((x1-x3)^2+(y1-y3)^2+(z1-z3)^2);
    c=sqrt((x3-x2)^2+(y3-y2)^2+(z3-z2)^2);
    S=Airetri(a,b,c);
    Vol=Vol+S*Rmoy/3;
end