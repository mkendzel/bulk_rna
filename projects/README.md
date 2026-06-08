## Project repository for all experiments

Each experiment gets its own folder here, containing that experiment's `data/`, `figures/`,
`results/`, and scripts (mirrors the layout used in `scRNA_pipeline`).

Experiments can be named for data sources (Plasmidsaurus codes or other company codes used for
alignment), for collaborations based on data generators, or for specific goals (funding/publication
opportunities).

All data in this repository must directly relate to bulk RNA analysis, alignment, and graph
generation.

### Setting up a new project

Copy `_template/` and rename it after your experiment:

```powershell
Copy-Item -Recurse projects\_template projects\YOUR_PROJECT_NAME
```

Shared resources (gene sets, gene-ID maps) live in the top-level `genesets/` folder and are
referenced from project scripts as `genesets/...`.
