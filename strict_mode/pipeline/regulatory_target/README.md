# Regulatory Targets

| Stage | Script | Purpose |
| --- | --- | --- |
| 06 | `06_define_regulatory_targets.py` | Intersect promoter-bound genes with DEGs and render regulatory-target figures. |
| 07 | `07_validate_outputs.py` | Verify figure inventory, DEG counts, and regulatory-target membership. |

Run this branch and its prerequisites through the parent workflow:

```bash
python ../run_pipeline.py --step regulatory --threads 16
```

The `figures_rendering/` directory contains required Python and R rendering
modules used by Stage 06. They are implementation code and are not run as
separate stages.
