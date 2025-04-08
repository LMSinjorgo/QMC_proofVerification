function [cGk_val] = computeEigenValueBounds(edgeIndicator,s)
% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
% compute c(G,s) = lambda_max(H_G)/2 - |E(G)| for all graphs in
% edgeIndicator, on s vertices

numGraphs = size(edgeIndicator,2);

% precompute all local hamiltonians so that we can quickly constract the
% hamiltonians corresponding to graphs.
% However, some edges never appear in the graphs, so we don't have to
% compute the local Hamiltonians H_e corresponding to every edge (minor
% speedup)
edgeOccurences = full(sum(edgeIndicator,2));
impossibleEdges = find(edgeOccurences == 0);
edgeIndicator(impossibleEdges,:) = [];
allEdges = uint16(nchoosek(1:s,2));
allEdges(impossibleEdges,:) = [];
[~,localHamiltonians] = createHamiltonians(allEdges);

% precompute values for the loop
cGk_val = zeros(1,numGraphs);
cardinalityEG = sum(edgeIndicator,1);
numPossibleEdges = uint16(size(allEdges,1));
oldEdgeIndicator = sparse(false(numPossibleEdges,1));
prevBound = 999;
H = sparse([],[],[],2^s,2^s,(2^s)*max(cardinalityEG));

% loop over all graphs and compute lambda_max(H_G)
for loopIdx = 1:numGraphs

    currentEdgeIndicator = edgeIndicator(:,loopIdx);
    diff = currentEdgeIndicator-oldEdgeIndicator;

    % We compute upper bounds on lambda_max(H_G), based on lambda_max(H_G)
    % of the previous graphs. If these upper bounds already prove c(G,k) <=
    % floor(s/2), we don't have to compute lambda_max(H_G)
    if nnz(diff) + prevBound <= floor(s/2)
        cGk_val(loopIdx) = nnz(diff) + prevBound;
        continue;
    else
        % check if the newly added edges form a star
        % if so, these edgesInStar can increase bound at most 1
        newEdges = allEdges(diff == 1,:);
        isStar = edgesInStar(newEdges);
        if isStar && ~isempty(newEdges)
            newBound = prevBound+1+nnz(diff == -1);
            if newBound <= floor(s/2)
                cGk_val(loopIdx) = newBound;
                continue;
            end
        else
            % if the new edges do not form a star graph, we can still
            % derive an upper bound as follows:
            % the new edges form a graph on 'numel(unique(newEdges(:)))'
            % vertices. The graph defined by these new edges has cGk at
            % most floor( numel(unique(newEdges(:))) / 2)
            % For each deleted edge, we also increase the bound by 1.
            newBound = prevBound + nnz(diff == -1) + floor(numel(unique(newEdges(:)))/2);
            if newBound <= floor(s/2)
                cGk_val(loopIdx) = newBound;
                continue;
            end
        end 
    end

    % if the previous bounds do not prove c(G,k) <= floor(s/2), we simply
    % construct H_G, and compute lambda_max(H_G). We can construct H_G
    % recursively based on H_G of the previous graph
    for e = 1:numPossibleEdges
        if diff(e) == 0
            continue;
        elseif diff(e) == 1
            H = H + localHamiltonians{e};
        elseif diff(e) == -1
            H = H - localHamiltonians{e};
        else
            error("programming error")
        end
    end

    % use the eigs function to compute the lambda_max(H_G)
    % and use lambda to compte c(G,k), with k = s.
    % Note that lambda_max(H_G) = eigs(H,1);
    prevBound = eigs(H,1)/2-cardinalityEG(loopIdx);
    cGk_val(loopIdx) = prevBound;

    % reset the old edge indicator
    oldEdgeIndicator = currentEdgeIndicator;
end
end

