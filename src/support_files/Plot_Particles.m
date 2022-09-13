function Plot_Particles(geode,Domain,Centres,Grains,VolumesGrains,FileName)

load(geode)
load(Domain)

figure1=figure;
set(figure1,'WindowStyle','docked')
for i=1:size(Centres,1)
    Centre=Centres(i,:);
    %if Centre(1,2)<0
    %if sqrt((Centre(1,1)-0.5)^2+(Centre(1,2)-0.5)^2+(Centre(1,3)-0.5)^2)>0.5
    %if -Centre(1,2)+Centre(1,3)-(1-Centre(1,1))>0
    %    continue
    %end
    Xgrain=real(Grains{i,1}(:,1));
    Ygrain=real(Grains{i,1}(:,2));
    Zgrain=real(Grains{i,1}(:,3));
    Rgrain=real(Grains{i,1}(:,4));
    VolumeGrain=real(VolumesGrains(i,1));
    %Color=Translate_Color(log(VolumeGrain),log(min(VolumesGrains)),log(max(VolumesGrains)));
    Color=Translate_Color(sqrt(VolumeGrain),sqrt(min(VolumesGrains)),sqrt(max(VolumesGrains)));
    %Color=Translate_Color(log(VolumeGrain),log(0.000005),log(0.0005));
    %Color=Translate_Color(Centre(1,2),0,1);
    %Color=rand(1,3);
    %surf1=trisurf(Structure,Xgrain,Ygrain,Zgrain,Zgrain,'LineStyle','none','BackFaceLighting','reverselit');shading interp;
    surf1=trisurf(Structure,Xgrain,Ygrain,Zgrain,'FaceColor',Color,'LineStyle','none','BackFaceLighting','reverselit');
    set(surf1,'AmbientStrength',.3,'DiffuseStrength',.3,'SpecularStrength',.1,'SpecularExponent',100);
    hold on
    %text(Centres(i,1),Centres(i,2),Centres(i,3),num2str(i),'FontSize',6)
end
Plot_Domain(Points_Domain,Struc_Domain)
xlabel('x (pixel)')
ylabel('y (pixel)')
zlabel('z (pixel)')
axis equal
grid off
camlight headlight
camlight right
camlight left
set(gcf,'Renderer','zbuffer')
%plot3(Centres(:,1),Centres(:,2),Centres(:,3),'*')
%figure1=figure;
%set(figure1,'WindowStyle','docked')
%plot(sort(VolumesGrains./VolumesCells),1/size(VolumesGrains,1):1/size(VolumesGrains,1):1)

set(gcf, 'InvertHardCopy', 'off');
File=[FileName,'.png'];
print('-dpng',File,'-r600')