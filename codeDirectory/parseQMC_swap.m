function [M,h,basis,cons] = parseQMC_swap(n)
% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
% Parse the QMC SDP relaxation in the swap algebra, corresponding to graphs
% on n vertices. The variables h_ij satisfy -s_ij, for the s the swap
% variable.
%
% It is known that the first level of the QMC SDP relaxation in the swap algebra is less tight than the second level of the QMC SDP 
% relaxation in the pauli matrices, see e.g., Section 5.1.2 in:
% 		Jun Takahashi, Chaithanya Rayudu, Cunlu Zhou, Robbie King, Kevin Thompson, and
% 		Ojas Parekh. An SU(2)-symmetric semidefinite programming hierarchy for quantum
% 		max cut. arXiv preprint arXiv:2307.15688v2, 2023. 
% 
% We take the algebraic relations of the s variables from
%
% Parekh, Ojas, and Kevin Thompson.
% "An optimal product-state approximation for 2-local quantum hamiltonians with positive terms."
% arXiv preprint arXiv:2206.08342 (2022).
% Section 4.3 (Monogamy of Entanglement on a Triangle)
%
% In particular, these relations state:
% (S_{ij})^2 = 1
% S_{ij} * S_{ij} = (1/2) * ( S_ij + S_ik + S_jk - 1)
% And these relations induce linear constraints on the moment matrix M
% (e.g., diag(M) == 1).
% Our basis used in the SDP relaxation is given by [1,h_12,..,h_1n,
% h_23,...,h_2n,....,h_{n-1,n}]

% Each index i in [n] is associated the i-th prime number, so that it is simple
% to uniquely identify the variables h_{ij} as p(i)*p(j), where p(i) is the
% i-th prime number.
primeList = primes(500);
primeList = uint32(primeList(1:n));

% create indices corresponding to the h-variables. In 'basis', each row
% corresponds to an edge, and the columns indicate the vertices. We uniquely identify
% h_{ij} with p(i)*p(j), for p(i) the i-th prime.
basis = uint32(nchoosek(1:n,2));
h_Primes = prod(primeList(basis),2,"native");
basisSize = numel(h_Primes);

% create the h variables in YALMIP.
h = sdpvar(basisSize,1);

% create the prime mapping, a vector that takes as input a prime,
% and returns the index of the corresponding h variable. For example,
% pimeMapping(6) = 1, because 6 = 2*3 = p(1)*p(2) and h_{12} is the first
% variable in the basis. Similarly, primeMapping(10) = 2, because 10 = 2*5
% = p(1)*p(3), and h_{13} is the second variable in the basis.
primeMapping = sparse(h_Primes,uint32(1),1:numel(h_Primes),max(h_Primes),1);

% create the moment matrix M in YALMIP, and add constraints that diagonal
% of M should equal 1. Also link the h-variables and M together, by
% constraining M(2:end,1) == h.
M = sdpvar(basisSize+1,basisSize+1);
cons = [M >= 0; M(2:end,1) == h; diag(M) == 1];

% We now compute the outer product of h_Primes (the primes that correspond
% to the h variables), to later construct equality constraints in M.
primeProducts = uint32(round(double([1;h_Primes])*double([1;h_Primes]')));

% For future use, create the matrix that maps 3 primes to the 3 unique products of 2 of the
% primes, i.e., [a,b,c] to [a*b,a*c,b*c]
mappingMat = [1,2;1,3;2,3];

% loop over upper triangular entries of the moment matrix M,
% and add additional equality constraints
for rowIdx = 2:size(M,1)-1
    for colIdx = rowIdx+1:size(M,2)


        % check the product of the primes/h variables.
        % If the variables correspond to disjoint edges, the product of
        % their corresponding primes should consist of 4 distinct prime
        % factors. If the h variables correspond to adjacent edges, there
        % should be 3 distinct prime factors.
        currentEntry = primeProducts(rowIdx,colIdx);
        primeFactors = factor(currentEntry);
        uniqueFactors = unique(primeFactors);
        numUniqueFactors = numel(uniqueFactors);

        if numUniqueFactors == 3
            % if the h variables share a vertex, this entry of the moment
            % matrix should be constrained via
            % S_{ij} * S_{ij} = (1/2) * ( S_ij + S_ik + S_jk - 1)
            %                 = (1/2) * (-h_ij - h_ik - h_jk - 1)
            % Thus, we need to determine the 3rd edges that creates the
            % triangle, for which we use the mappingMat.
            allEdges = uniqueFactors(mappingMat);
            edgeProducts = prod(allEdges,2,'native');

            % now h_Indices contains the linear indices of the edges in the
            % triangle
            h_Indices = full(primeMapping(edgeProducts));

            cons = [cons; M(rowIdx,colIdx) == (-1-sum(h(h_Indices)))/2];
        end

    end
end
end