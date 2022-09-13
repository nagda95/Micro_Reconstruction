function [PD,E,V] = pwrDiagramPD(T, PC)
% function PD = pwrDiagramPD(T, PC)
%
% T: triangulation
% PC: power centers of triangles in T
%
% The output cell PD contains pieces of the power diagram, indexed by
% dimension. PD{mP} contains pieces of dimension zero (power centers that
% correspond to fully-dimensional pieces of the triangulation T) and PD{1}
% contains fully-dimensional regions of the power diagram (corresponding to
% vertices of the triangulation).

P = piecesPD(T);
mP = size(P,1);
PD = cell(mP,1);

% Do I need this? Yes! I use PD{4}
for i=[2,4]
    EA = edgeAttPD(T, P{i}); 
    for j=1:size(EA,1)
        EA{j} = PC(EA{j},:);
    end
    PD{i} = EA;
end
clear EA

i=1;
V = edgeAttPD(T, P{i}); % vertex indices for each cell
for j=1:size(V,1)
    EA{j} = PC(V{j},:);
end
PD{i} = EA;
clear EA

i=3;
E = edgeAttPD(T, P{i}); % vertex indices for each edge
for j=1:size(E,1)
    EA{j} = PC(E{j},:);
end
PD{i} = EA;