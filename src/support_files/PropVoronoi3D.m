function [Volumes,Voisins]=PropVoronoi3D(X,Y,Z,Vertices,Cellules)
Volumes=zeros(size(Cellules,1),1);
Voisins=zeros(size(Cellules,1),1);
for n=1:size(Cellules,1)
    C=transpose(Cellules{n,1});
    P=Vertices(C,:);
    K=convhull(P+(10^-10)*randn(size(P,1),3));
    [Pd,Sd]=Tri2Domain(P,K);
    [Volume]=Volume_Domain(Pd,Sd);
    Volumes(n,1)=Volume;
end
Voisins=Voisinage3D(X,Y,Z,delaunay(X,Y,Z));