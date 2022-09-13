function [PD, PDinf,P,C,E,Neigh,infCells,Edges] = powerDiagramWrapper_mod(X, wts,scaling)
% This is a modification of powerDiagramWrapper by Frederick McCollum 
% (https://www.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams)
% Modifications by David Bourne, Heriot-Watt University, 4 July 2021
%
% function [PD, PDinf] = powerDiagramWrapper(X, wts)
%
% X: set of points
% wts: weights of the points in X
%
% The output cell PD contains the pieces of the power diagram indexed by
% dimension. PD{1} contains fully-dimensional regions of the power diagram.
% PDinf contains points on infinite edges of the power diagram.
%
% If the initial points are in R^2, the power diagram is drawn.

% PD{1} contains the coordinates of the vertices of the NONEMPTY cells
% PD{2} contains the coordiantes of the vertices of the faces
% PD{3} contains the coordinates of the vertices of the edges
% PD{4} contains the coordinates of the vertices
% infCells contains the indices of the unbounded cells
% Neigh{i} contains the indices of the neighbours of cell i

%disp('Lifting points');
LX = liftPD(X, wts);
%disp('Computing convex hull');
K = convhulln(LX,{'Qt'});
%disp('Extracting lower hull');
ind = normalsPD(LX, K);
T = K(ind, :);
clear K LX ind

nT = size(T,2);
%disp('Finding pieces');
[P, total] = piecesPD(T);  % Why compute total ????????????????
nonEmpty = P{1}; % indices of non-empty cells
%disp('Finding power centers');
[PC, powers] = powercentersPD(T, X, wts);  % Why compute powers ????????????
%disp('Finding free boundary');
FF = freeBouPD(T, P{nT-1});
%disp('Finding power diagram');
[PD,E,V] = pwrDiagramPD(T, PC); % Why recompute P inside this function? Pass it as an argument ??????????
C = cell(size(X,1),1);
C(nonEmpty)=V;
clear V

%V = PD{4};

% Do I need this?
infCells=unique(FF)'; % list of unbounded cells %%%%%%%%%%%%%%%%%%%%%%%%%

% Can I delete the following?
Neigh=cell(size(X,1),1); % store the neighbour relations, Neigh{j} = indices of neighbours of cell j, %%%%%%%%%%%%%%%%%%%%
for j=1:size(P{2},1) % for each face in the diagram
     Neigh{P{2}(j,1)}(end+1)=P{2}(j,2);
     Neigh{P{2}(j,2)}(end+1)=P{2}(j,1);
end

% edge indices for each cell
Edges=cell(size(X,1),1); % store the edges in each cell, Edges{j} = indices of edges in cell j %%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:size(P{3},1) % for each edge in the diagram
    Edges{P{3}(j,1)}(end+1)=j;
    Edges{P{3}(j,2)}(end+1)=j;
    Edges{P{3}(j,3)}(end+1)=j;
end

% % Makes list of vertices: Verts{j} = list of vertices in cell j
% % The ordering of the vertices is the same as PD{4}, i.e., 
% % PD{4}(Verts{i})=cell containing the coordinates of the vertices in cell i
% Verts = cell(size(X,1),1);  
% numVerts = size(P{4},1);
% for j=1:numVerts
%     Verts{P{4}(j,1)}(end+1)=j;
%     Verts{P{4}(j,2)}(end+1)=j;
%     Verts{P{4}(j,3)}(end+1)=j;
%     Verts{P{4}(j,4)}(end+1)=j;
% end

center = mean(X,1);
PDinf = zeros(size(FF));

% find distance from center to farthest powercenter
%length = max(sqrt(sum(bsxfun(@minus, PC, center).^2,2)));
length = max(scaling,max(sqrt(sum(bsxfun(@minus, PC, center).^2,2)))); %%%%%%%%%%%%%%%

%disp('Finding points on infinite edges of the power diagram');
for i=1:size(FF,1)
    facet = X(FF(i,:),:);
    ea = edgeAttPD(T, FF(i,:));
    pc = PC(ea{1},:);
    ct = mean(facet,1);
    
    % find vector normal to the facet
    v = null(bsxfun(@minus, facet(1,:), facet(2:end,:)))';
    
    % reorient v to point outward
    if dot(center - ct, v) > 0
        v = -1*v;
    end
    
    % scale v to ensure newpt is sufficiently far away
    v = length*v;
    
    % find point on infinite edge of power diagram
    newpt = pc + v;
    p = piecesPD(FF(i,:));
    for j=1:size(p,1)
        for k=1:size(p{j},1)
            ind = find(ismember(P{j},p{j}(k,:), 'rows'));
            PD{j}{ind} = [PD{j}{ind}; newpt];
        end
    end
    
    % keep track of generated point
    PDinf(i,:) = newpt;
end