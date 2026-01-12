\### Format for Stabilizer Check Matrix / Normalizer Generator Matrix



We adopt MacKay’s nonbinary alist format to represent both stabilizer

check matrices and normalizer generator matrices.



Each matrix corresponds to an M × N operator matrix, where:

\- N is the number of physical qubits (columns),

\- M is the number of generators (rows).



---------------------------------------

Pauli Encoding

---------------------------------------

Each nonzero matrix entry is encoded as an integer representing the Pauli

operator acting on a given qubit:



0 : Identity (I)

1 : X

2 : Z

3 : Y



(All unspecified positions are implicitly I.)



---------------------------------------

Alist Structure

---------------------------------------

The file is organized as follows:



1\. N M

&nbsp;  (number of columns, number of rows)



2\. max\_col\_weight max\_row\_weight



3\. Column weights (length N)

&nbsp;  wt(c1), wt(c2), ..., wt(cN)



4\. Row weights (length M)

&nbsp;  wt(r1), wt(r2), ..., wt(rM)



5\. Column descriptions (N lines)

&nbsp;  For each column c:

&nbsp;    (row\_index, Pauli\_type) pairs for all non-identity entries in column c,

&nbsp;    padded with zeros to length max\_col\_weight.



6\. Row descriptions (M lines)

&nbsp;  For each row r:

&nbsp;    (column\_index, Pauli\_type) pairs for all non-identity entries in row r,

&nbsp;    padded with zeros to length max\_row\_weight.



---------------------------------------

Example: \[\[5,1,3]] Stabilizer Code

---------------------------------------



Stabilizer check matrix:



X Z Z X I

I X Z Z X

X I X Z Z

Z X I X Z



Corresponding alist representation:



5 4

4 4

3 3 3 4 3

4 4 4 4



(col descriptions)

1 1 3 1 4 2 0 0

1 2 2 1 4 1 0 0

1 2 2 2 3 1 0 0

1 1 2 2 3 2 4 1

2 1 3 2 4 2 0 0



(row descriptions)

1 1 2 2 3 2 4 1

2 1 3 2 4 2 5 1

1 1 3 1 4 2 5 2

1 2 2 1 4 1 5 2









