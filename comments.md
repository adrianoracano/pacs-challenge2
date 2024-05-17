Comments in the source indicated with @note

Nice the idea of a main for the comparison. But to compare the result is is different of making the norm of the difference, its not necessary showing all elements!

The makefile deandeas not correspond to that indicated in the README: you have to chenge it to compile the different sources! You could have just set differrnt targets! I have modified
the makefile to show how different targets can be handles. I have also colelcted all object files that do not contain a main() in a static library. This is a good practice to avoid recompiling.

