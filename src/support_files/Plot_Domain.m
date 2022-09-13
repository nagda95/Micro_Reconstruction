function Plot_Domain(Points_Domain,Struc_Domain)
%figure
for i=1:size(Struc_Domain,1)
    xv=Points_Domain(Struc_Domain{i,1}(:,1),1);
    yv=Points_Domain(Struc_Domain{i,1}(:,1),2);
    zv=Points_Domain(Struc_Domain{i,1}(:,1),3);
    plot3(xv,yv,zv,'-b')
    hold on
end
axis equal