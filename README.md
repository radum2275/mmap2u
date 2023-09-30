# mmap2u
Credal Marginal MAP solver version 1.0

# Building the project
This is a standard `cmake` C++ project so please feel free to use your favourite build toolchain. It should work out of the box with `vs code` and the `cmake` extensions.

# Running the solver
Once the binary is compiled, the solver can be run with the following command line arguments. Note that the default name of the compiled binary is `loopy` (you can change it in the `CMakeLists.txt` file).

The benchmark problems are assumed to be given in the UAI format.

```
./loopy --input-file alarm.uai --query-file alarm.query --task MMAP --algorithm hc --query-type maximax --seed 12345678 --iterations 10 --threshold 0.0000001 --flip-proba 0.2 --alpha 0.9 --ibound 2 --max-flips 10000 --taboo-size 100 --cache-size 1000000 --time-limit 3600
```

The command line arguments are:

* --input-file : the problem file in UAI format (e.g., problem.uai)
* --query-file : the query file (e.g., problem.query)
* --algorithm : the CMMAP algorithm (cve - CVE, dfs - DFS, hc - SHC, ts - TS, sa - SA)
* --query-type : the CMMAP task (maximax, maximin)
* --iterations : the number of iterations used by SHC, TS, SA
* --max-flips : maximum number of flips used by SHC, TS, SA
* --flip-proba: the flip probability used by SHC
* --alpha : the cooling schedule used by SA
* --taboo-size : the size of the taboo list used by TS
* --ibound : the ibound used by CMBE(i)
* --time-limit : the time limit in seconds (-1 means no limit)