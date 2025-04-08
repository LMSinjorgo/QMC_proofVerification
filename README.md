
# Computer proof of lemma $c(G,k) \leq \lfloor s/2 \rfloor$

## By S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
A computer verification of Lemma 3.10 in the paper:

*Sander Gribling, Lennart Sinjorgo, and Renata Sotirov. 
Improved approximation ratios for the Quantum Max-Cut problem on general, triangle-free and bipartite graphs.* 

## The verification procedure
Here, we provide computational details on the verification of $c(G,s) \leq \lfloor s /2 \rfloor$ for all $G \in \mathcal{G}_{s}$ for $s \in \{5,7,9,11,13\}$. Afterwards, we provide more details on how to run said code and the dependencies. 

We performed the following steps on a laptop (16 GB RAM and Intel i7-1165G7 CPU), which required approximately 16 hours to run. 

We first use the software package $\texttt{nauty}$ [[2]](#2) to generate the graphs in $\mathcal{G}_s$ that satisfy the properties 1 to 3 of Lemma 3.9 in the paper. Then, for the case $s \in \{5,7,9\}$, we verify that $c(G,2) \leq \lfloor s /2 \rfloor$ for these graphs in $\mathcal{G}_s$ using SDP. Computing $c(G,2)$ can be done using the Pauli-based SDP relaxation, but, in fact, solving an SDP relaxation of  based on the SWAP operators already suffices to establish the desired upper bound. (More precisely, we compute the first level of the QMC SDP relaxation based on the SWAP operators, see  [[1,3]](#1) and in particular the discussion in  [[3, Sect. 5.1.2]](#3).   

For $s = 11$, we consider the 26360 (see Table 3.2) remaining graphs $G_j = ([11],E_j)$, $j \in [26360]$ in the sequence as returned by $\texttt{nauty}$. We construct the corresponding Hamiltonians recursively as

$$H_{G_{j+1}} = H_{G_{j}} + \sum_{e \in E_{j+1} \setminus E_j} H_e - \sum_{e \in E_{j} \setminus E_{j+1}} H_e$$

Constructing $H_{G_{j+1}}$ using this recursion is efficient since $\texttt{nauty}$ returns a sequence of graphs where $E_j \approx E_{j+1}$. Additionally, this recursion shows that 

$$    \lambda_\mathrm{max}\left(H_{G_{j+1}}\right) \leq \lambda_\mathrm{max}\left(H_{G_{j}}\right) + \lambda_\mathrm{max}\left( \sum_{e \in E_{j+1} \setminus E_{j}} H_e \right).$$

Here, we have used that $H_e = (1/4)H_e^2 \succeq 0$. Generally, $\lambda_\mathrm{max}\left( \sum_{e \in E_{j+1} \setminus E_j} H_e \right) \leq 4 \left| E_{j+1} \setminus E_j \right|$, and tighter bounds are possible if, for example, the edges $E_{j+1} \setminus E_j$ form a star graph. If the above bound already proves that $c \left( G_{j+1},s \right) = \lambda_\mathrm{max}  \left( H_{G_{j+1}}\right)/2 - | E_{j+1} | \leq \lfloor s /2 \rfloor$, we do not carry out the computation of $\lambda_\mathrm{max} \left( H_{G_{j+1}}\right)$, nor the construction of $H_{G_{j+1}}$.  If the above bound does not prove $c(G_{j+1},s) \leq \lfloor s /2 \rfloor$, we compute $\lambda_\mathrm{max}\left( H_{G_{j+1}}\right)$ with the $\texttt{eigs}$ function in MATLAB. 

The case $s = 13$ proceeds similarly as the case $s = 11$, except we first discard some of the 9035088 graphs in 
$G_{13}$ that satisfy properties 1 to 3 of Lemma 3.9, by arguing as follows: of these 9035088 graphs, 959842 of them do not satisfy property 4 with $|S| = 2$. We verify that these 959842 graphs do not satisfy property 4 in approximately 5 seconds. Of the now remaining 8075246 graphs, 1622184 of them satisfy $\tau(G) \leq 6$. We verify that these 1622184 graphs satisfy $\tau(G) \leq 6$ in approximately 20 seconds. By Item 2 of Lemma 3.8, these graphs satisfy $c(G,13) \leq \tau(G) \leq 6$, so we may discard them. For the remaining 6453062 graphs $G$, we compute upper bounds on $\lambda_\mathrm{max} \left(H_G \right)$ in the manner described for the case $s = 11$.

## Dependencies
- [Gurobi](https://www.gurobi.com/)
- [Mosek](https://www.mosek.com/)
- [nauty](https://pallini.di.uniroma1.it/)
- [YALMIP](https://yalmip.github.io/)

This code uses Gurobi for solving ILPs and MOSEK for solving SDPs. It is also possible to use alternative ILP/SDP solvers, although the code will have to be slightly adjusted.

## Usage of code
(see the files in codeDirectory)

Run the $\texttt{nauty}$ commands from $\texttt{nautyCommands.txt}$, to obtain the files humanReadable_s[X].txt, where $X \in \{ 7,9,11,13\}$ (the case $X=5$ does not require $\texttt{nauty}$). These files contain the graphs $G \in  \mathcal{G}_s$ for which we need to verify that $c(G,k) \leq  \lfloor s/2  \rfloor$.

Then run verifyLemma\_s5\_to\_s13.mat. This code performs the full verification (which requires approximately 16 hours of computation time). Most of this time is spent on the case $s = 13$. The cases $s \in \{5,7,9,11\}$ require approximately 5 minutes. The codes verify\_s[X].mat, for $X \in \{5,7,9,11,13\}$, verify the cases of $X$ separately.

## References
<a id="1">[1]</a>  Adam Bene Watts, Anirban Chowdhury, Aidan Epperly, J. William Helton, and Igor
Klep. Relaxations and exact solutions to quantum max cut via the algebraic structure
of swap operators. Quantum, 8:1352, 2024.

<a id="2">[2]</a>  Brendan D. McKay and Adolfo Piperno. Practical graph isomorphism, II. Journal of
Symbolic Computation, 60:94–112, 2014.

<a id="3">[3]</a>  Jun Takahashi, Chaithanya Rayudu, Cunlu Zhou, Robbie King, Kevin Thompson, and
Ojas Parekh. An SU(2)-symmetric semidefinite programming hierarchy for quantum
max cut. arXiv preprint arXiv:2307.15688v2, 2023.
