## Hydropower Dams Pareto Optimization

Code for running Pareto Optimization for Hydropower dam planning projects.

Based on the original work in [1][2], executes a dynamic programming recursive algorithm
for finding the exact or approximate Pareto frontier for hydropower dam planning across
multiple objectives. Updated and optimized for [3][4], which include algorithmic
optimizations [3] and compression and expansion for running subsets of objectives
in parallel, as well as new methods for identifying solution representations from [TODO: Delia paper]
. See the respective papers for algorithmic details.

### Library Requirements

This version of the code requires boost and OpenMP libraries.

See https://www.boost.org/doc/libs/1_84_0/more/getting_started/ for installing Boost.

If compiling on Mac, you may need to download the openmp libraries separately, see https://mac.r-project.org/openmp/

### Compiling

To compile the code, run:

```
cmake CMakeLists.txt
make
```

which will compile an executable program HydroDams in the same directory as the code

### Input File Format

If using a non-provided input file, it must follow the following format:

```
p #num_nodes #num_edges #num_criteria
criteria crit_1 crit_2 ... crit_k
d dam_1_id dam_1_status dam_1_crit_val_1 dam_1_crit_val_2 ... dam_1_crit_val_k
...
d dam_m_id dam_m_status dam_m_crit_val_1 dam_m_crit_val_2 ... dam_m_crit_val_k
n node_1_id node_1_crit_val_1 node_1_crit_val_2 ... node_1_crit_val_k
...
n node_n_id node_n_crit_val_1 node_n_crit_val_2 ... node_n_crit_val_k
r node_id
e edge_1_parent_node_id edge_1_child_node_id edge_1_dam_id
...
e edge_m_parent_node_id edge_m_child_node_id edge_m_dam_id
```

where:
* num_nodes should equal num_edges + 1
* the criteria list has k = #num_criteria names
* each of the dam and node lines contains d objective values in the same order
as the criteria information
* there should be m = #num_edges lines for the dams (lines starting with *d*)
* there should be n = #num_nodes lines for the nodes (lines starting with *n*)
* dam_status should be the the status of the dam, 0 or 1 (not built or built)
* dam_crit_val is a triple (s|p|q) where is the value when building, p is the
passage probability for building, and q is the passage probability for not building
* node_crit_val is just the value associated with the node
* the *r* line designates the root of the tree (mouth of the river)
* the *e* lines are adjacency lists that connect parent to child nodes via a specific dam,
that is it defines how the dam splits up the river to upstream and downstream sections

### Running

After compiling the code, to run a simple test, simply execute

```
./HydroDams
```

This will run a baseline optimization for the entire Amazon with
just two objectives: energy and connectivity. The executable takes a number of configurable
parameters. To see the help information, run

```
./HydroDams -h
```

Default parameters, such as sorting method and affine transform inclusion, are set
to the optimal parameters based on [3].

Important parameters to set include:

- -epsilon: epsilon approximation value, it is not recommended to use 0 beyond 3 objectives
- -path: the path to the data file, if not using the provided files

### Solutions

Once a run is complete, a solution.txt file is created in the results directory (by default
this will be results/run_name/solution.txt). The format of this file starts
with a series of lines containing the meta data for the given run, ending with
"criteria considered" (the list of objectives being optimized) and "criteria all"
(the list of all criteria calculated). All lines after this is the list of solutions in the Pareto frontier
with a comma separated list first of all the criteria values in the same order as the
"criteria all" information, then the number of dams built in that solution, followed
by a space separated list listing out all the dam ids that were built for that solution.

### Compression and Expansion

Based on work in [4], for especially large problems, such as running 6+ objectives for the
entire Amazon river basin, it is infeasible to calculate the Pareto frontier for
reasonable epsilon values. As such, for these large problems it is recommended to
instead run in parallel a number of subsets of objectives and to combine the resulting
Pareto frontiers into a single frontier. This has weaker theoretical guarantees, but is
much more efficient to calculate.

The basic principle is to take linear combinations of different objectives to produce
a smaller number of objectives that are optimized against. Thus you convert *d* objectives to *k*
objectives, where each of the *k* objectives is a linear combination of one or more
of the *d* objectives, mutually exclusive from each other. To run a single instance with a
subset of objectives, run:

```
./HydroDams -lp -w k k_1 <weights> ... k_k <weights>
```

where *k* is the number of objectives to compress to, *k_1* is the number of objectives
in the first compressed objective, followed by a list of linear weights within that objective,
which should all sum to 1. The order the objectives are compressed is the order they are listed
in for the -criteria parameter.

**Energy objectives should not be grouped with non-energy objectives.**

For example, say we have the 5 objectives: energy, connectivity, biodiversity, ghg, and dor.
We can compress these 5 objectives into 3 objectives, leaving energy by itself,
and grouping connectivity with biodiversity and ghg with dor as follows:

```
./HydroDams -criteria 5 energy connectivity biodiversity ghg dor -lp -w 3 1 1.0 2 0.5 0.5 2 0.5 0.5
```

Note in this example, within their groups, each objective is weighed equally, but
to get good coverage of the Pareto frontier, it is recommended you run different
parallel runs with a mixture of different weights, keeping in mind the weights must add
up to 1.0 within each group.


### Papers

[1] X. Wu et al., “Efficiently Approximating the Pareto Frontier:
Hydropower Dam Placement in the Amazon Basin,” AAAI, vol. 32, no. 1, Apr. 2018, doi: 10.1609/aaai.v32i1.11347.

[2] J. M. Gomes-Selman, Q. Shi, Y. Xue, R. García-Villacorta, A. S. Flecker, and C. P. Gomes,
“Boosting Efficiency for Computing the Pareto Frontier on Tree Structured Networks,” in Integration of
Constraint Programming, Artificial Intelligence, and Operations Research, vol. 10848, W.-J. van Hoeve, Ed.,
in Lecture Notes in Computer Science, vol. 10848. , Cham: Springer International Publishing, 2018, pp. 263–279.
doi: 10.1007/978-3-319-93031-2_19.

[3] \<COMING SOON\> M. Grimson et al., "Scaling Up Pareto Optimization for Tree Structures with Affine Transformations:
Evaluating Hybrid Floating Solar-Hydropower Systems in the Amazon," AAAI, 2024

[4] Y. Bai, Q. Shi, M. Grimson, A. Flecker, and C. P. Gomes, “Efficiently Approximating High-Dimensional Pareto
Frontiers for Tree-Structured Networks Using Expansion and Compression,” in Integration of Constraint Programming,
 Artificial Intelligence, and Operations Research, vol. 13884, A. A. Cire, Ed., in Lecture Notes in Computer Science,
  vol. 13884. , Cham: Springer Nature Switzerland, 2023, pp. 1–17. doi: 10.1007/978-3-031-33271-5_1.