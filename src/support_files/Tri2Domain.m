function [Points_Domain,Struc_Domain]=Tri2Domain(Points,Structure);
Points_Domain=Points;
Struc_Domain=cell(size(Structure,1),1);
for i=1:size(Structure,1)
    Struc_Domain{i,1}=Structure(i,:)';
end