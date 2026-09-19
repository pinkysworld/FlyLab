# Architecture

```
compound + concentration
        |
        v
  occupancy engine     (Hill / Langmuir, cited EC50)
        |
        +-----------------------------+
        v                             v
 insect receptor map          vertebrate panel
 (nAChR, RDL, GluCl)         (nAChR subtypes, GABA-A)
        |
        v
  weight / gain patch on named synapses of a connectome subgraph
        |
        v
  circuit runtime (LIF or rate)
        |
        v
  assay readout (MN9 Hz, DNp01 Hz, …) + notebook JSON
```

The connectome is a **netlist**. Pharmacology is a **patch** on that netlist. Keep them separate files so a reviewer can see which number came from Janelia and which number came from a paper EC50.
