function [isStar] = edgesInStar(edgeList)
if isempty(edgeList)
    isStar = false;
    return
end
% check if all edges in edgelist form a star graph
numEdges = size(edgeList,1);
if numEdges == 1
    isStar = true;
    return;
end

% let the first edge be {v1,v2}. Then check if all edges are adjacent to
% v1, or if all edges are adjacent to v2.
for vertexCounter = 1:2
    v = edgeList(1,vertexCounter);
    % first check v1
    for k = 2:numEdges
        newEdge = edgeList(k,:);
        if any(newEdge == v)
            if k == numEdges
                % all edges are adjacent to v
                isStar = true;
                return;
            end
            continue;
        else
            break;
        end
    end
end

isStar = false;
end

