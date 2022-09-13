function [V,C] = croppedLaguerre3D(x,y,z,w,L1,L2,L3,c)
% Create a 3D Laguerre diagram with seeds with coordinates (x,y,z) and
% weights w, cropped to a rectangular box with sides of length L1, L2, L3
% and centre c.
%
% This function relies on a modification of powerDiagramWrapper. The
% original function was written by
% Frederick McCollum (https://www.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams)
%
% Author: David Bourne, Heriot-Watt University, 4 July 2021

N = length(x); % number of generators

% make sure x,y,z,w,c are column vectors
x = reshape(x,N,1);
y = reshape(y,N,1);
z = reshape(z,N,1);
w = reshape(w,N,1);
c = reshape(c,3,1);
X=[x,y,z];

% Create the unbounded Laguerre diagram
[PD,PDinf,P,C,E,Neigh,infCells,Edges] = powerDiagramWrapper_mod([x,y,z],w,2*(norm(c)+norm([L1,L2,L3])));

% Crop to the rectangular box

% Add each corner of the box to the correct cell: For each corner x, just compute its power_i(x) for each cell i and take the min.
corners = repmat(c,1,8)+[[L1;L2;L3],[L1;L2;-L3],[L1;-L2;L3],[L1;-L2;-L3],[-L1;L2;L3],[-L1;L2;-L3],[-L1;-L2;L3],[-L1;-L2;-L3]]/2;
for i=1:8 % for each corner
    powers = zeros(N,1);
    for k=1:N % for each cell %%% Could be speeded up by only taking k such that the cell does not lie entirely outside the box %%%%%%%%% - but shouldn't happen in Tata problem
       powers(k)=dot(corners(:,i)'-X(k,:),corners(:,i)'-X(k,:))-w(k);
    end
    [~,cellInd] = min(powers); % find which Laguerre cell the corner belong to 
    % Add the corner as a vertex of the cell
    cornerInd = length(PD{4})+1; % index of the new vertex
    PD{4}{end+1}=corners(:,i)'; % add the coords on the new vertex to the list of vertex coords
    C{cellInd}=[C{cellInd},cornerInd];
end

numE = length(E); % number of edges in the Laguerre diagram (including infinite edges)
edgeInd = 1:numE; % indices of edges

% linear constraints corresponding to the sides of the box: If x satisfies
% A*x-b < 0, then x lies strictly inside the box.
A = [eye(3,3);-eye(3,3)];
b = [c;-c]+[L1;L2;L3;L1;L2;L3]/2;

% intersect each edge of the power diagram with the bounding box

for i=1:6 % for each plane of the bounding box
    
    for j=edgeInd % for each edge of the Laguerre diagram
        
        v = PD{3}{j}'; % coordinates of the vectices at the end of the edge (note: one vertex may be a point on an infinite edge), dimension 3-by-2
        
        if(A(i,:)*v-b(i)>0) % if both ends of the edge lie outside the bounding plane
            
            edgeInd = setdiff(edgeInd,j); % delete the edge from the list of edges (not strictly necessary but saves some time)
            vertInd = E{j}; % get the indices of the vertices
            for k=1:length(vertInd) % length(vertInd)= 2 (usually) or 1 (if one of the vertices lies on an infinite edge)
                cellInd = P{4}(vertInd(k),:); % list of cells containing this vertex
                for m=cellInd % delete this vertex from every cell it belongs to
                    C{m}=setdiff(C{m},vertInd(k));
                end
            end
            
        elseif(A(i,:)*v-b(i)<0) % else if both ends of the edge lie inside the bounding plane
            
            % do nothing (note that the edge could still lie outside the
            % bounding box)
            
        else % else one end of the edge lies inside the bounding plane, and one end lies outside
            
            % delete the outside vertex from any cell in which is appears
            
            vertInd = E{j}; % get the indices of the endpoint vertices
            
            %%%%% ALSO NEED TO TREAT THE CASE WHERE THERE IS NO FINITE VERTEX?
            %%%%% And where the bounding box lies strictly inside one cell?
            if(length(vertInd)==1) % if there is only one finite vertex (i.e., if the edge is unbounded)
                
                vFinite = PD{4}{vertInd}'; % coordinates of the finite vertex
                if(A(i,1)*vFinite-b(i)>0) % if the finite vertex lies outside the plane:
                    cellInd = P{4}(vertInd,:); % list of cells containing this vertex
                    for m=cellInd % delete this vertex from every cell it belongs to
                        C{m}=setdiff(C{m},vertInd);
                    end
                end
                
            else % else the edge has two finite vertices
                
                % find which one is inside and which one is outside
                verts=PD{4}{vertInd(1)}';
                if(A(i,:)*verts-b(i)>0)
                    vOut=vertInd(1); % vertex index of the outside vertex
                    vIn=vertInd(2); % vertex index of the inside vertex
                else
                    vOut=vertInd(2); % vertex index of the outside vertex
                    vIn=vertInd(1); % vertex index of the inside vertex
                end
                cellInd = P{4}(vOut,:); % list of cells containing the outside vertex
                for m=cellInd % delete the outside vertex from every cell it belongs to
                    C{m}=setdiff(C{m},vOut);
                end
                
            end
            
            p = intersectLinePlane(v,A(i,:)',b(i)); % intersect the edge with the bounding plane
            
            if(isOnBoundary(p,A,b)) % if the intersection point lies on the boundary of the box, then
                
                % add the intersection point as a vertex to the correct cells
                
                pInd = length(PD{4})+1; % index of the new vertex
                PD{4}{end+1}=p'; % add the coords on the new vertex to the list of vertex coords
                cellInd = P{3}(j,:); % cells corresponding to edge j
                for m=cellInd % add the new vertex to the relevant cells
                    C{m}=[C{m},pInd];
                end
                
            else
                
                %                 % delelete the inner vertex if it is a finite vertex
                %                 if(length(vertInd)~=1)
                %                     cellInd = P{4}(vIn,:); % list of cells containing the inside vertex
                %                     for m=cellInd % delete the outside vertex from every cell it belongs to
                %                         C{m}=setdiff(C{m},vIn);
                %                     end
                %                 end
                %                 % delete the edge
                %                 edgeInd = setdiff(edgeInd,j);
                
                % delelete the inner vertex if it is a finite vertex
%                 if(length(vertInd)~=1)
%                     cellInd = P{4}(vIn,:); % list of cells containing the inside vertex
%                     for m=cellInd % delete the inside vertex from every cell it belongs to
%                         C{m}=setdiff(C{m},vIn);
%                     end
%                 else
%                     vFinite = PD{4}{vertInd}'; % coordinates of the finite vertex
%                     if(A(i,1)*vFinite-b(i)<0) % if the finite vertex lies inside the plane:
%                         cellInd = P{4}(vertInd,:); % list of cells containing this vertex
%                         for m=cellInd % delete this vertex from every cell it belongs to
%                             C{m}=setdiff(C{m},vertInd);
%                         end
%                     end
%                 end
                % delete the edge
                %edgeInd = setdiff(edgeInd,j);
                
            end
            
        end
        
    end
    
end

% intersect edge of the bounding box with each Laguerre cell. Only need to
% consider Laguerre cells that are not strictly inside the bounding box

edgesBB = {repmat(c,1,2)+[[-L1;-L2;L3],[-L1;L2;L3]]/2,...
    repmat(c,1,2)+[[L1;-L2;L3],[L1;L2;L3]]/2,...
    repmat(c,1,2)+[[L1;-L2;L3],[-L1;-L2;L3]]/2,...
    repmat(c,1,2)+[[L1;L2;L3],[-L1;L2;L3]]/2,...
    repmat(c,1,2)+[[-L1;-L2;-L3],[-L1;L2;-L3]]/2,...
    repmat(c,1,2)+[[L1;-L2;-L3],[L1;L2;-L3]]/2,...
    repmat(c,1,2)+[[L1;-L2;-L3],[-L1;-L2;-L3]]/2,...
    repmat(c,1,2)+[[L1;L2;-L3],[-L1;L2;-L3]]/2,...
    repmat(c,1,2)+[[L1;-L2;-L3],[L1;-L2;L3]]/2,...
    repmat(c,1,2)+[[-L1;-L2;-L3],[-L1;-L2;L3]]/2,...
    repmat(c,1,2)+[[L1;L2;-L3],[L1;L2;L3]]/2,...
    repmat(c,1,2)+[[-L1;L2;-L3],[-L1;L2;L3]]/2};

tol=1e-12; % to avoid rounding errors

for i=1:12 % for each edge of the bounding box
    
    v1=edgesBB{i}(:,1); % one end of the edge
    v2=edgesBB{i}(:,2); % the other end of the edge
    
    for k=1:N % for each nonempty cell  % it would be more efficient to do this for each face rather than each cell
        
        if(~isempty(C{k}))
            
            numVerts = length(C{k});
            verts = cell2mat(PD{4}(C{k}))';
            constraints = (A*verts-repmat(b,1,numVerts) < tol);
            
            %if(~constraints) % If the cell does not lie strictly inside the box (i.e., it touches the boundary) we might need to add more vertices
                
                An = zeros(length(Neigh{k}),3);
                bn = zeros(length(Neigh{k}),1);
                count=0;
                for j=Neigh{k} % build the contraints corresponding to the boundary of the cell
                    % Add a linear constraint for each neighbour (An*z <= bn)
                    count=count+1;
                    An(count,:) = X(j,:)-X(k,:);
                    bn(count) = (dot(X(j,:),X(j,:))-dot(X(k,:),X(k,:))+w(k)-w(j))/2;
                end
                
                % intersect the edge of the box with each face of the cell:
                count=0;
                for j=Neigh{k} % for each face of the cell
                    
                    count=count+1;
                    a = An(count,:);
                    b =bn(count);
                    
                    % check if the edge has one vertex either side of the face
                    if((a*v1-b)*(a*v2-b)<0)
                        
                        % intersect the edge and the face
                        p = intersectLinePlane([v1,v2],a,b);
                        
                        % check whether the intersection point p lies on the
                        % boundary of the cell
                        if(isOnBoundary(p,An,bn))
                            
                            % add the intersection point to the list of
                            % vertices of the cell
                            pInd = length(PD{4})+1; % index of the new vertex
                            PD{4}{end+1}=p'; % add the coords on the new vertex to the list of vertex coords
                            C{k}=[C{k},pInd];
                            
                        end
                        
                    end
                    
                end
                
            %end
            
        end
        
    end
    
end

V = cell2mat(PD{4});

end
 

function p = intersectLinePlane(l,a,b)
% Intersects a line and a plane.
% l contains the coordinates of the end points of the line, dimension 3-by-2.
% a, b specify the plane: The plane is given by the equation dot(a,x)-b = 0, 
% where a is a column vector (dimension 3-by-1) and b is a scalar.
    
v1 = l(:,1);
v2 = l(:,2);
t = (b-dot(a,v1))/dot(a,v2-v1);
p = v1 + t*(v2-v1);

end


function [out]=isOnBoundary(p,A,b)
% check whether point p lies on the boundary (or inside) of the box specified by A,b
tol = 1e-12;
if(A*p-b < tol)
    out=1;
else
    out=0;
end

end