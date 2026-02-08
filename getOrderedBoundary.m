function boundaryNodes = getOrderedBoundary(E)
% getBoundary - retrieves boundary of triangular element mesh
%   boundaryNodes = getOrderedBoundary(E)
%   retrieves a list of all nodes on the boundary of the element mesh,
%   ordered by adjacency
% Ben Sarfati 2/2026

%If tris are quadratic, order nodes clockwise as opposed to conventionally
[numElems,nodesPerElem] = size(E);
if nodesPerElem == 6
    E = E(:,[1 4 2 5 3 6]);
end

%Create list of all edges of all elements, remove duplicates (counting
%oppositely ordered edges as the same edge), and track which edges had
%duplicates
[edgesNoDupes,~,indOld2IndNew] = unique(sort([E(:) [E(numElems+1:end) E(1:numElems)]'],2),'rows');

%Create boundary edge list (edges that had no duplicates)
boundaryEdges = edgesNoDupes(accumarray(indOld2IndNew,1)==1,:);
numUniqueEdges = length(boundaryEdges);

%Create list of all nodes belonging to these edges ordered by adjacency
boundaryNodes = [boundaryEdges(1); zeros(numUniqueEdges-1,1)];
row = 1;
col = 1;
for i = 2:numUniqueEdges%+1 for last node = first node
    boundaryNodes(i) = boundaryEdges(row,3-col);
    boundaryEdges(row,:) = [];
    [row,col] = find(boundaryEdges == boundaryNodes(i),1);
end