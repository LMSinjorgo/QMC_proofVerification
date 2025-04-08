function [satisfiesStableSetCondition] = stableSetCondition(edgeIndicator)
% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
% stable set condition
% check for which graphs |N(S)| <= |S|. While this condition holds for
% general stable sets S, it is most efficient to restrict |S| = 2.
% (which is also what this code does)
% In that case, we look for two non-adjacent vertices of degree two, that are
% adjacent to the same two vertices. I.e., we look for a 4-cycle subgraph


[maxNumEdges,numGraphs] = size(edgeIndicator);
satisfiesStableSetCondition = false(1,numGraphs);
edgeIndicatorTranspose = edgeIndicator';

% compute s, the number of vertices, from the number of maxNumEdges
s = round(0.5*(1+sqrt(1+8*maxNumEdges)));

% compute the degrees of the vertices per vertex in each graph, thus given
% a (s times numGraphs) vertex.
[degrees] = computeDegrees(edgeIndicator,s);

% transpose for faster access
degrees = degrees';

% create the linear indices for edges
lowerTriuIdx = find(tril(ones(s),-1));
edgeMat = zeros(s); edgeMat(lowerTriuIdx) = 1:numel(lowerTriuIdx);
edgeMat = edgeMat+edgeMat';
% now edge ij will have linear index edgeMat(i,j)

fullEdges = nchoosek(1:s,2);
graphIndexOriginal = uint32(1:numGraphs);
for edgeCounter = 1:size(fullEdges,1)


    e = fullEdges(edgeCounter,:);
    v1 = e(1); v2 = e(2);
    edgeLinIdx = edgeMat(v1,v2);

    % select all graphs that have deg(v1) = deg(v2) = 2, ...
    graphSubsetBool = and(degrees(:,e(1)) == 2,degrees(:,e(2)) == 2);
    % ...and {v1,v2} not an edge. The linear index of {v1,v2} is given by
    % edgeLinIdx.
    graphSubsetBool = and(graphSubsetBool,~edgeIndicatorTranspose(:,edgeLinIdx));

    % select the graphs that satisfy the necessary properties
    graphSubset = edgeIndicator(:,graphSubsetBool);
    
    % get the numbers/indices of the graphs in the subset
    graphIndex = graphIndexOriginal(graphSubsetBool);

    % check if vertices are common. To do so, first get the indicent edges
    % (the linear index)
    [neighbors1,neighbors2] = getNeighborEdges(v1,v2,edgeMat);

    graphSums = and(graphSubset(neighbors1,:),graphSubset(neighbors2,:));
    commonVertices = sum(graphSums,1);

    isSeparable = commonVertices == 2;
    originalIndex = graphIndex(isSeparable);

    satisfiesStableSetCondition(originalIndex) = true;
end

end

function [degrees] = computeDegrees(edgeIndicator,s)
numGraphs = size(edgeIndicator,2);
degrees = zeros(s,numGraphs);

allEdges = nchoosek(1:s,2);
for vertex = 1:s
    adjacentEdges = (allEdges(:,1) == vertex) + (allEdges(:,2) == vertex);
    adjacentEdges = adjacentEdges > 0;

    degrees(vertex,:) = sum(edgeIndicator(adjacentEdges,:),1);
end
end

function [neighbors1,neighbors2] = getNeighborEdges(v1,v2,edgeMatrix)
% get the edge indices of edges adjacent to either one of the vertices
% in edge. Note that the vertices in V = {v1, v2} are not connected.
% neighborsi returns the vertices adjacent to vi
edgeMatrix(v1,v2 ) = 0;
edgeMatrix(v2,v1 ) = 0;

neighbors = edgeMatrix(:,[v1,v2]);
neighbors1 = neighbors(:,1);
neighbors2 = neighbors(:,2);
neighbors1(neighbors1 == 0) = [];
neighbors2(neighbors2 == 0) = [];
end




