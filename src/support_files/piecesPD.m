function [P, total] = piecesPD(T)
% function [P, total] = piecesPD(T)
%
% T: triangulation
%
% The output P contains all pieces of the triangulation, from vertices
% (dimension zero) to fully-dimensional facets (dimension n-1).

%%%% In 3D:
%%%% P{1} <-> nonempty cells of the power diagram: P{1} = list of indices of nonempty cells.
%%%% P{2} <-> faces of the power diagram: P{2}(m,:) = [i,j] means that i and j are neighbours and share a face
%%%% P{3} <-> edges of the power diagram: P{3}(m,:) = [i,j,k] means that cells i,j,k meet along the edge m
%%%% P{4} <-> vertices of the power diagram: P{4}(m,:) = [i,j,k,l] means that cells i,j,k,l meet at the vertex m

[m, n] = size(T);
P = cell(n,1);
total = 0;

P{1} = unique(T);
P{1} = P{1}(:);
total = total + size(P{1},1);

for i=2:n-1
    Q = [];
    for j=1:m
        %Q = [Q; combnk(T(j,:),i)];
        Q = [Q; nchoosek(T(j,:),i)];
    end
    Q = sort(Q, 2);
    P{i} = unique(Q, 'rows');
    total = total + size(P{i},1);
    clear Q;
end

P{n} = T;
total = total + size(P{n},1);