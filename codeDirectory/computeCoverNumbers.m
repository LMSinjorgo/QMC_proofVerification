function [coverNumsUB] = computeCoverNumbers(edgeIndicator)
% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
% given the edgeIndicator matrix (graphs correspond to columns, and edges
% to rows), check whether the bound obtained by the cover number already
% proves c(G,k) <= floor(s/2).
%
% Note that tau(G) can be computed by an Integer linear programming problem
% (ILP). Instead of solving this ILP to compute tau(G), we simply check
% whether the ILP defined by tau(G) <= floor(s/2) is feasible. (because we
% don't care about the precise value of tau(G), only whether tau(G) <=
% floor(s/2) )

[maxNumEdges,numGraphs] = size(edgeIndicator);

% compute s, the number of vertices, from the number of maxNumEdges
s = round(0.5*(1+sqrt(1+8*maxNumEdges)));
allEdges = uint16(nchoosek(1:s,2));
orderedEdgeList = allEdges'; orderedEdgeList = orderedEdgeList(:);

% construct the matrix A. The polytope defined by Ax <= [1;...,1; floor(s/2)] defines the
% feasibility problem of tau(K_s) <= floor(s/2), where K_s is the complete
% graph on s vertices. We will select the rows of A that correspond to the
% edges in specific graphs later.
% The bottom row of A is all ones, and corresponds to the constraint
% sum(x) <= floor(s/2)
A = sparse(uint16(repelem(1:maxNumEdges,2)),orderedEdgeList,1,maxNumEdges,s);
A = [A; ones(1,s)];


coverNumsUB = zeros(1,numGraphs,"uint16");
numEdges = full(sum(edgeIndicator,1));
edgeIndicatorPlusTrue = [edgeIndicator; true(1,numGraphs)];

% set gurobi parameters
model.vtype = repmat('B',s,1);
params.outputflag = 0;
geqSense = repmat('>',maxNumEdges,1);
oldCover = zeros(s,1);

% test if gurobi is installed
try
    model.A = A ;
    model.sense = [geqSense(1:size(A,1)-1); '<'];
    model.rhs = [ones(size(A,1)-1,1);floor(s/2)];
    result = gurobi(model, params);
catch gurobiError
    error("Gurobi not (properly) installed!")
end

for loopIdx = 1:numGraphs
    currentEdgeIndicator = edgeIndicatorPlusTrue(:,loopIdx);
    currentNumEdges = numEdges(loopIdx);
    newA = A(currentEdgeIndicator,:);

    % try the cover previously found and see if it is also a cover for this
    % graph
    if all(newA*oldCover >= 1)
        coverNumsUB(loopIdx) = floor(s/2);
        continue;
    end

    model.A = newA;
    model.sense = [geqSense(1:currentNumEdges); '<'];
    model.rhs = [ones(currentNumEdges,1);floor(s/2)];

    result = gurobi(model, params);
    if result.mipgap == 0
        % if the ILP is feasible, tau(G) <= floor(s/2)
        coverNumsUB(loopIdx) = uint16(floor(s/2));
        oldCover = result.x > 0.5;
        if nnz(oldCover) > floor(s/2)
            error("Cover is too large!")
        end
    else
        % if it is not feasible, we use the upper bound tau(G) <= s
        coverNumsUB(loopIdx) = uint16(s);
    end


end

end

