# DC-opf model
We modified the original 'qps_ot.m' file to use the built-in 'dual-simplex' algorithm as a fallback when the 'interior-point' algorithm fails to solve the linear optimization problem. 
**PLEASE REPLACE THE 'qps_ot.m' FILE IN MATPOWER v7.1 WITH THE VERSION PROVIDED IN THIS FOLDER.**
Keep in mind that using different solvers may produce varying results.

## Precautions
The built-in solver in MATPOWER v7.1 sometimes encounters difficulties locating feasible solutions, particularly within the IEEE300 benchmark. 
To address this issue, we relax the minimal power generation requirement when switching on the generator, allowing $P_i^{+}$ to change continuously from 0 to its maximum value. 
Consequently, the mixed integer optimization is transformed into a linear optimization problem that is often easier to solve. 
Moreover, we ignore any load that cannot be dispatched by MATPOWER v7.1. 
The resulting optimal (or minimal) total power loss $l_p^{(\min)}$ is rounded to 8 decimal places.

