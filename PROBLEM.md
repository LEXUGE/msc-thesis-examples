We work over GF(2) and consider binary matrices.

# Goal

Given an arbitrary binary **target matrix with linearly independent rows**

$$
T \in \mathbb F_2^{r\times n},
$$

find a matrix $M=[A\mid B]$ satisfying the constraints below.  The target
code is its row space $C_T=\operatorname{rowspan}(T)$. Only this row space
matters: replacing $T$ by a row-equivalent matrix does not change the
problem. The supplied rows themselves must be linearly independent; redundant
or zero target rows are not valid input.

The number $n$ of B-columns and the number $r$ of target rows are arbitrary;

## The matrix

The parameters $k$ and $m$ are free. The matrix

$$
M \in \mathbb F_2^{k\times(m+n)}
$$

is written as $M=[A\mid B]$, where

$$
A\in\mathbb F_2^{k\times m},\qquad B\in\mathbb F_2^{k\times n}.
$$

Every row of $M$ has Hamming weight 3: exactly two `1`s occur in its A-part
and exactly one `1` occurs in its B-part. Every A-column must have positive
weight; equivalently, the graph represented by A has no isolated vertices.

## Constraint 0

There must be a basis $b_1,b_2,\ldots$ of $\ker A$ whose supports partition
the A-columns. In other words, the support of each basis vector is disjoint from
the others, and together the supports contain every A-column.

For example, if $m=5$, then

$$
b_1=(1,1,0,0,0),\quad b_2=(0,0,1,1,0),\quad b_3=(0,0,0,0,1)
$$

satisfies this condition.

## Constraint 1

Consider every row combination $z^TM$ whose A-part vanishes:

$$
z^TA=0.
$$

The corresponding B-parts must span exactly the target code:

$$
\{z^TB : z^TA=0\}=C_T.
$$

Also, no nonzero zero-A row combination may have a zero B-part:

$$
z^TA=0\ \text{and}\ z^TB=0 \quad\Longrightarrow\quad z=0.
$$

Equivalently, the map from zero-A row combinations to their B-parts is
injective. Since the rows of $T$ are independent, this gives

$$
\dim\ker(A^T)=r.
$$

Equivalently, after eliminating the A-part by GF(2) row operations, the
remaining B-rows may be any basis of the same row space as $T$. They do not
need to equal the supplied rows of $T$ verbatim.

## Constraint 2

The basis from Constraint 0 partitions the A-columns into groups. For every
group and every weight-2 vector $v\in\mathbb F_2^m$ supported in that group,
there must be a $w\in\mathbb F_2^k$ such that

$$
w^TM=[v\mid v']
$$

for some $v'\in\mathbb F_2^n$ with Hamming weight at most 2.

## Optimisation objective

For a fixed target matrix $T$, find a feasible $M$ while minimizing $k$,
the number of rows of $M$. The value of $m$ may vary as part of the search.

The intended next step is an automated search, for example with OR-Tools, that
accepts an arbitrary target matrix $T$ and either produces such an $M$ or
proves infeasibility within specified search bounds.

## Concrete example targets

You are provided with two example targets.

n = 16, r = 11
```
[[0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0]
 [0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0]
 [0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1]
 [0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 1]
 [0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1]
 [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1]
 [0 1 1 0 0 0 0 0 1 0 0 1 0 0 0 0]
 [1 0 1 0 0 0 0 0 0 0 0 0 1 0 1 0]
 [0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0]
 [1 0 1 0 0 0 0 0 0 1 0 1 0 0 0 0]
 [0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 1]]
```

n = 18, r = 10
```
[[1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0]
 [0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1]
 [0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0]
 [0 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0]
 [0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 0 0]
 [0 1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0]
 [0 0 0 0 0 0 1 1 0 0 1 0 1 0 0 1 0 0]
 [0 0 1 0 1 0 0 0 0 0 1 0 0 1 0 0 0 1]
 [1 0 1 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0]
 [1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 1 0]]
```

As a first step you should try find solutions to the two concrete examples above, make sure to understand what helped you solve the problem, what is generalizable etc.
