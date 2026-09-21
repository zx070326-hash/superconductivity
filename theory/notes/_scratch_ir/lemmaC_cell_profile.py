import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from lemmaC_cell_partition import dk_ratio_t, trial_at
cells=[(0.09,0.36,4),(0.36,1.44,6),(1.44,4.32,16),(4.32,8.64,12),(8.64,17.28,16),(17.28,34.56,24),(34.56,69.12,48),(69.12,72.25,24)]
for xa,xb,Ns in cells:
    t=trial_at(np.sqrt(xa*xb),Ns)
    xs=np.exp(np.linspace(np.log(xa),np.log(xb),17))
    r=[dk_ratio_t(t,x) for x in xs]
    print(f"[{xa:7.3f},{xb:7.3f}] Ns={Ns:2d}  ratio: " + " ".join(f"{v:.2f}" for v in r))
