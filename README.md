# InAsSb_atomic_structures
To investigate the effect of Sb incorporation into the InAs crystal structure 
at the atomic level, we used artificial neural network based potentials.


https://github.com/inholeegithub/InAsSb_atomic_structures.git

## Workflow Overview

```mermaid
flowchart TD
    A[Build InAs Unit Cell] --> B[Introduce Sb Atoms]
    B --> C[Load Pretrained ANN Potential]
    C --> D[Fully Relax Unit Cell]
    D --> E[Run Atomic-Scale Simulations]
    E --> F[Analyze Sb Incorporation Effects]
```
