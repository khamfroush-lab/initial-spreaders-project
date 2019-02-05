%%%
%%% Last modified by Injung Kim, 3/23/2015
%
% Small Workd Network (Watts & Strogatz Model)
% A graph with n*deg/2 edges is constructed, i.e. the nodal degree is n*deg for
% every node
%
%  INPUTS:
%   n: number of nodes of the graph to be generated, n>=1
%   deg: mean degree, deg%2=0 && 0<deg<n-1
%   rp: rewiring probability, 0<=rp<=1
%  Output:
%   adj - adjacency matrix representing the generated graph
%  
function adj = small_world_graph(n, deg, rp)
% Construct a regular lattice: a graph with n nodes, each connected to deg
% neighbors, deg/2 on each side.
dHalf = deg/2;
rows = reshape(repmat([1:n]', 1, deg), n*deg, 1);
columns = rows+reshape(repmat([[1:dHalf] [n-dHalf:n-1]], n, 1), n*deg, 1);
columns = mod(columns-1, n) + 1;
B = sparse(rows, columns, ones(n*deg, 1));
A = sparse([], [], [], n, n);

% With probability rp rewire an edge avoiding loops and link duplication.
% Until step i, only the columns 1:i are generated making implicit use of A's
% symmetry.
for i = [1:n]
    % The i-th column is stored full for fast access inside the following loop.
    col= [full(A(i, 1:i-1))'; full(B(i:end, i))];
    for j = i+find(col(i+1:end))'
        if (rand()<rp)
            col(j)=0;
            deg = randi(n);
            while deg==i || col(deg)==1
                deg = randi(n);
            end
            col(deg) = 1;
        end
    end
    A(:,i) = col;
end

% A is not yet symmetric: to speed things up, an edge connecting i and j, i < j
% implies A(i,j)==1, A(i,j) might be zero.
T = triu(A);
A = T+T';

% initiate adj matrix, then put the scale free network
adj = zeros(n);
adj = full(A);

end 