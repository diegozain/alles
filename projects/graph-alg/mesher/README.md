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
```
---

### Example

[![](../pics/mesh.png)](./)

