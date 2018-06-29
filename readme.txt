GF2mat.T_fact(T,U,P) is fast for full rank binary matrix, but relatively slow for matrix with smaller rank. The function defined in itpp is rather short. I may found or build some other program to run it.(reference?)


threshold of Bacon Shor code.

Pick size 9x9=81,25x25=625,51x51=2601,125x125=15625

error probability p=(0.010,0.050,0.002)

number of error: n_e=1,000,000

fraction of failures: n_bad/n_e

note: the previous calculation going very slow. To get quick result with use of multiple CPU, I should run the task seperately.

steps: fromg G get S first and save S into a file. Then use S to calculate with different p in different process.


12/04/2017:
next step
1. show shreshold plot for bacon shor code
2. do same calculation on generalized bacon shor code and plot the shreshold
3. combine belief propagation method and random window decoder for toric codes.

group meeting:
check by hand: input error: single error
syndrome: only two neighboring generators unsatisfied
the bp decoding should converge correctly in this case.

1/31/2018:
start using tmux for everything.
run bp decoding for toric codes with 3 different size. get the plot of P_succ versus p. The default value of Xi is zero. Then try some negative value like -10.-20. The way to optimize Xi is to find the minimum average weight of the out put error from the BP decoding in the unconveged vased.
If I need to see the type of error that doesn't convege for BP decoding, I need to print them in the actual pattern (square lattice)


3/10/2018:
make a Z-error-only decoder. change the difinition in code_generate.c, remove the tilde part in decoding. 

updated files:
code_generate.c for toric code
bp_decoding3.c

3/12/2018
moved eoutdated code files into a sub folder

