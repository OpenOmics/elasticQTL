# elasticQTL

> A reproducible QTL pipeline for PLINK-based association testing, LD clumping, nested elastic net modeling, stability filtering, and final variant annotation — plus an external scoring/refit workflow for WGS/WES/array cohorts.

---
## Overview

`elasticQTL` provides two main workflows:

- **Training**: run association testing, LD clumping, genotype export, nested elastic net modeling, stability filtering, and final variant annotation
- **Validation / external scoring**: apply a trained model to external cohorts and optionally refit on the shared marker set before scoring

Use the entrypoints:

- `pipeline` for training
- `validate` for external scoring/refit

For full details, see:


- [`docs/PIPELINE.md`](docs/PIPELINE.md) — training workflow steps, nested CV, missing genotype policies, and outputs
- [`docs/VALIDATION.md`](docs/VALIDATION.md) — external cohort scoring, refit workflow, and cohort harmonization'
- [`docs/INPUT_FORMATS.md`](docs/INPUT_FORMATS.md) — required inputs, config fields, and manifest formats
- [`docs/METHODS.md`](docs/METHODS.md) — workflow details, modeling approach, and missing genotype handling
- [`docs/TROUBLESHOOTING.md`](docs/TROUBLESHOOTING.md) — common errors and fixes


  
## Quickstart

### Managed install (Skyline / OpenOmics)

The pipeline is installed at:

```
/data/openomics/prod/elasticQTL/v0.1.0/
```

> ⚠️ **Do not run on the head node.** The pipeline will error immediately if attempted.

```bash
# Grab an interactive node
srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash

# Add executables to PATH
export PATH="/data/openomics/prod/elasticQTL/v0.1.0/bin:${PATH}"

# View help menus
pipeline -h
validate -h
```

If run on the head node, you will see:

```
############################################################
ERROR: Running this pipeline on a head node is not allowed.
Please submit this as a job to the cluster or run it from
an interactive node. You can grab an interactive node with:

  srun -N 1 -n 1 --time=8:00:00 --mem=64gb -c 4 --pty bash
############################################################
```

---
### Local clone

If you are setting up elasticQTL outside of a managed install, clone the repo directly:

```bash
# Clone
git clone https://github.com/your-org/elasticQTL.git
cd elasticQTL

# Make entrypoints executable
chmod +x bin/pipeline bin/validate

# Add to PATH (add this to your ~/.bashrc or ~/.bash_profile to persist)
export PATH="$(pwd)/bin:${PATH}"
```

Then verify the install:

```bash
pipeline -h
validate -h
```
---

### Running the pipeline
```bash
# 1. Copy and edit a training config
cp config/params.template.env config/study.env
nano config/study.env

# 2. Run training
pipeline --config config/study.env

# 3. Copy and edit a validation config
cp config/validation.template.env config/validation.env
nano config/validation.env

# 4. Run external scoring/refit
validate --config config/validation.env
```


## Contributing

PRs and issues are welcome. When filing a bug report, please include:

- The config used *(redact paths if needed)*
- The relevant `logs/*.log` files
- `manifest/params_used.txt` (Step 1) or the validation logs (Step 2)
