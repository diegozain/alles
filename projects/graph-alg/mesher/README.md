# Mesh on a grid
diego domenzain

April 2021 @ Colorado School of Mines

## Mesh with neighbor information

__Given a rectangular gridded mesh where only a certain region is wanted for modeling, how do we extract this wanted region?__

This code answers this question by finding three key constructs:

1. find only the nodes that matter in the domain, 
1. find their neighbors,
1. find what type of neighbors they are.

The words **mesh** or **mesh-grid** refer to the existing rectangular gridded mesh.

The word **graph** refers to the region of interest inside the mesh-grid.

The *Matlab* code uses lots of ```for, while``` and ```if``` **on purpose**. It is meant to be **pseudocode** for *Fortran*.

```
1. The only nodes that matter in the domain are given in these two column vectors,

graph2mesh : indexes are graph nodes, entries are mesh nodes
mesh2graph : indexes are mesh nodes, entries are graph nodes

graph2mesh:

| |
| |
| | size = # of nodes that matter by 1
| |
| |
| |

mesh2graph:

| |
| |
| | size = # of nodes in the mesh by 1
| |
| |
| |
| |
| |

2. Their neighbors are given by these two matrices,

neigh_mesh  : row indexes are graph nodes. Row entries are neighbors of that node, in the mesh.
neigh_graph : row indexes are graph nodes. Row entries are neighbors of that node, in the graph.

3. Their neighbor type is given by this matrix,

neigh_type : row indexes are graph nodes. Row entries are the type of neighbor for that node.

neigh_mesh, neigh_graph, and neigh_type are all of size:

|        |
|        |
|        | size = # of nodes that matter by max # of neighbors in the mesh
|        |
|        |
|        |

since we are assuming a 2D mesh-grid, max # of neighbors in the mesh is 4.

the way the neigh_mesh, neigh_graph, and neigh_type refer to the neighbors of node 'i' (i.e. in row 'i') is by,

     2
     |
3 -- i -- 1
     |
     4

that is, columns 1, 2, 3, and 4 represent neighbors right, up, left and down.
```

The last steps are to build the PDE operator ```L```.

```
we need L to be defined as a sparse matrix,
so we need arrays I and J to do:

       L = sparse(I,J,v);

where 'v' is special.

1. we need the total number of non-zero entries on L.
that is for each row i of L, we have as many entries as neighbors of i + 1.
the +1 is to count its own entry.

2. we also need to count for each node in the graph, 
how many entries in I (and J) are of that node.

for 1, we sum all positive entries of neigh_type + n_g2m.
the term n_g2m accounts for each L(i,i) entry.

for 2, we sum columnwise all positive entries of neigh_type.

n_IJ : total number of non-zero entries of L (also length of I and J).
n_ij : holds the info of how many entries in I (and J) belong to each node i.
       it is an array of size n_g2m by 1.

now we need to build I and J.

both I and J are of size n_IJ by 1.

I denotes the row entries, and J the column entries.

for each node i:
   I needs to have n_ij(i)+1 consecutive entries with the number i.
   In the same place as these n_ij(i)+1 consecutive entries, 
      J needs to have the number i in the first entry, 
      and then each neighbor of i (in the graph) in the subsequent entries.
```
---

### Example

[![](../pics/mesh.png)](./)

