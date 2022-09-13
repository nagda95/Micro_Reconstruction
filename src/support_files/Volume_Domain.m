function [Volume]=Volume_Domain(Points_Domain,Struc_Domain)
Volume=0;
XD=mean(Points_Domain(:,1));
YD=mean(Points_Domain(:,2));
ZD=mean(Points_Domain(:,3));
for iface=1:size(Struc_Domain,1)
    XA=Points_Domain(Struc_Domain{iface,1}(1,1),1);
    YA=Points_Domain(Struc_Domain{iface,1}(1,1),2);
    ZA=Points_Domain(Struc_Domain{iface,1}(1,1),3);
    if size(Struc_Domain{iface,1},1)>3
        for itri=2:size(Struc_Domain{iface,1},1)-2
            XB=Points_Domain(Struc_Domain{iface,1}(itri,1),1);
            YB=Points_Domain(Struc_Domain{iface,1}(itri,1),2);
            ZB=Points_Domain(Struc_Domain{iface,1}(itri,1),3);
            XC=Points_Domain(Struc_Domain{iface,1}(itri+1,1),1);
            YC=Points_Domain(Struc_Domain{iface,1}(itri+1,1),2);
            ZC=Points_Domain(Struc_Domain{iface,1}(itri+1,1),3);
            Volume=Volume+Tetravol(XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD);
        end
    elseif size(Struc_Domain{iface,1},1)==3
        XB=Points_Domain(Struc_Domain{iface,1}(2,1),1);
        YB=Points_Domain(Struc_Domain{iface,1}(2,1),2);
        ZB=Points_Domain(Struc_Domain{iface,1}(2,1),3);
        XC=Points_Domain(Struc_Domain{iface,1}(3,1),1);
        YC=Points_Domain(Struc_Domain{iface,1}(3,1),2);
        ZC=Points_Domain(Struc_Domain{iface,1}(3,1),3);
        Volume=Volume+Tetravol(XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD);
    end
end