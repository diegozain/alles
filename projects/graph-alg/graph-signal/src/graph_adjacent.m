function V = graph_adjacent(ngraph,n_ij,n_IJ,I,J,W)
% diego domenzain
% december 2021
% ------------------------------------------------------------------------------
%
% discretization of
%         the laplacian of a graph 🍇 encoded in I & J & W
%
% where
%       L = sparse(I,J,V);
% ------------------------------------------------------------------------------
% I & J & W
%
% both I and J are of size n_IJ by 1.
% I denotes the row entries, and J the column entries.
%
% for each node i:
% I needs to have n_ij(i)+1 consecutive entries with the number i.
% in the same place as these n_ij(i)+1 consecutive entries,
% J needs to have the number i in the first entry,
% and then each neighbor of i (in the graph) in the subsequent entries.
% ------------------------------------------------------------------------------
%
% the ith row of matrix L has the following form,
%
%    L(i_g2m,:) = [   0        -w(of edge i⏤j) ]
%                  ith entry       jth entries
%
% ------------------------------------------------------------------------------
% ngraph     • total # of nodes in the graph.
% n_ij       • holds the info of how many entries in I (and J) belong to each
%              node i. It is an array of size ngraph by 1.
% n_IJ       • total number of non-zero entries of L (also length of I and J).
% I          • I denotes the row entries
% J          • J the column entries
% ------------------------------------------------------------------------------
%    L(i_g2m,:) = [   0        -w(of edge i⏤j) ]
%                  ith entry       jth entries
% ------------------------------------------------------------------------------
V = zeros(n_IJ,1);
% ------------------------------------------------------------------------------
% begining of row in L
il = 1;
% end of row in L (initialized)
il_= 0;
% ------------------------------------------------------------------------------
% ⚫ access a node 'i_g2m' in the graph,
%    which is indexed by 'il' in I, J, & V.
for i_g2m=1:ngraph
    % -- end of row in L
    il_ = il_ + n_ij(i_g2m)+1;
    % -- in row of L
    % ◼ loop thru neighbors of 'i_g2m',
    %    whose position in I, J, & V is indexed by 'ij',
    %    → asign values to,
    %
    %                    L(i_g2m , ij) = - weight(of edge i_g2m⏤ij)
    %
    %    → and do the neighbor sum for the entry L(i_g2m,i_g2m),
    %
    %                        L(i_g2m,i_g2m) = Σj weight(of edge i_g2m⏤ij)
    %    inside this loop:
    %      • J(il) gives the ith node in the graph,
    %      • J(ij) gives neighbor of ith node in the graph.
    % --
    for ij=(il+1):il_
        % entry of L(i_g2m , ij)
        V(ij) = W(ij);
    end
    % -- begining of next row in L
    il = il_ + 1;
end
end
