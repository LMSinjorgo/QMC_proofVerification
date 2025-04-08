function [H,localHamiltonianCell] = createHamiltonians(edgeList,W)
% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
% input: edge List (i.e., [3,4 ; 2,5] if graph has edges {3,4} and {2,5}
%           edge weight matrix W
% output: matrix H_G = sum_{e \in E} H_e, and the local Hamiltonians H_e, according to the QMC problem.
% if edgelist is empty return zero
if isempty(edgeList)
    H = 0; localHamiltonianCell = []; return;
end

numVertices = max(max(edgeList));

[I,x,y,z,S] = createPauli();
x = sparse(x);
z = sparse(z);
y = sparse(y);
I = sparse(I);

if nargin == 1
    W = ones(size(edgeList,1),1);
end

H = sparse([],[],[],2^numVertices,2^numVertices);
cellInput = cell(1,numVertices);
localHamiltonianCell = cell(1,size(edgeList,1));
for edgeIdx = 1:size(edgeList,1)
    
    % create the kron Product
    [cellInput{:}] = deal(I);

    currentEdge = edgeList(edgeIdx,:);
    
    cellInput(currentEdge) = {x,x};
    localHamiltonian = speye(2^numVertices)-kronProd(cellInput);

    cellInput(currentEdge) = {y,y};
    localHamiltonian = localHamiltonian-kronProd(cellInput);

    cellInput(currentEdge) = {z,z};
    localHamiltonian = localHamiltonian-kronProd(cellInput);

    localHamiltonian = localHamiltonian;
    localHamiltonianCell{edgeIdx} =W(edgeIdx)* localHamiltonian;
    H = H + W(edgeIdx) * localHamiltonian;   
end




end

function [A] = kronProd(cellInput)
numCellElements = numel(cellInput);

A = kron(cellInput{end-1},cellInput{end});
if numCellElements == 2
    return;
end


for p = numCellElements-2:-1:1
    currentMat = cellInput{p};
    A = kron(currentMat,A);
end
end

function [I,x,y,z,S,H] = createPauli()
% output the pauli matrices
x = ones(2)-eye(2);
y = [0, -1i;
    1i,0];
z= diag([1,-1]);

I = eye(2);
S = eye(4);
S(2:3,2:3) = [0,1; 1,0];
H = kron(I,I) - kron(x,x) - kron(y,y) - kron(z,z);
end

