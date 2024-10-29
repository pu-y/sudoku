int: order;
int: nrow = order * order;
int: ncol = order * order;
set of int: rows = 1..nrow;
set of int: cols = 1..ncol;
array[rows, cols] of int: grid;

array[rows, cols] of var 0..nrow: fill;
array[rows, cols] of var 1..nrow: completed;

constraint forall(i in rows, j in cols)(completed[i, j] = fill[i, j] + grid[i, j]);
constraint forall(i in rows, j in cols)(grid[i, j] != 0 -> fill[i, j] = 0);

include "alldifferent.mzn";
constraint forall(i in rows)(alldifferent(completed[i, ..]));
constraint forall(j in cols)(alldifferent(completed[.., j]));
constraint forall(i in 0..((nrow div order) - 1), j in 0..((ncol div order) - 1))
                (alldifferent([completed[order*i + u, order*j + v] | u in 1..order, v in 1..order]));

solve satisfy;

output ["\(completed[i, j])" ++ if ((j mod ncol) = 0) then "\n" else " " endif  |i in rows, j in cols];