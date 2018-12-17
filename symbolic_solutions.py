## Hui's ORC expected richness
from sympy import *
R, C1, c2, c3 = symbols ( 'R C1 c2 c3')
f1 = C1 + c2*R + c3*log(R)
solve(f1, R)

## Csi value from observed and total species in Tovo's TNB
S, p, k, C = symbols ( 'S p k C')
f1 = (1-(1-p)**k) * C - S
solve(f1, p)
