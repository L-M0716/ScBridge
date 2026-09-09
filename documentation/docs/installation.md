<div class="installation-guide" markdown="1">

# Installation

!!! info "Installation overview"
    ScBridge packages its analysis software in an Apptainer container.
    Installation consists of obtaining the workflow, building the image, and verifying that it runs.

## Requirements

| Requirement | Notes |
|---|---|
| Apptainer | Required on a compatible execution host to run the ScBridge SIF image |
| ScBridge repository | Provides the launcher, workflow rules, configuration, and build definition |
| Disk space | Allow space for input data, references, intermediate objects, and results |

R, Python, Snakemake, STAR (including STARsolo), and analysis packages are supplied through the
container rather than installed individually on the host. The commands below use a Bash-compatible
shell on the machine running Apptainer.

## 1. Check Apptainer

```bash
apptainer --version
```

If the command is unavailable, install Apptainer following its
[official installation guide](https://apptainer.org/docs/admin/latest/installation.html),
or use the runtime provided by your computing facility. A container image does not replace
the host container runtime.

## 2. Get the workflow

Download the ScBridge source code and enter the repository directory:

```bash
git clone https://github.com/L-M0716/ScBridge.git
cd ScBridge
```

## 3. Build the image

Run from the repository root, where `run_workflow`, `workflow1/`, and `build/` are located:

```bash
apptainer build --fakeroot scRNA_seq.sif \
  build/scRNA_seq.def
```

!!! note "Build requirements"
    The build needs network access to retrieve the base image and dependencies, sufficient space
    for temporary files and the final SIF image, and support for `--fakeroot` on the execution host.
    Keep the repository files referenced by the definition available during the build.

## 4. Verify the runtime

```bash
apptainer exec scRNA_seq.sif echo "Container started successfully"
apptainer exec scRNA_seq.sif snakemake --version
apptainer exec scRNA_seq.sif STAR --version
apptainer exec scRNA_seq.sif R --version

apptainer exec \
  -B "$PWD":/opt/scRNA_workflow \
  scRNA_seq.sif \
  bash /opt/scRNA_workflow/run_workflow --help
```

Run these commands from the repository root. `$PWD` is the host repository directory, mounted at
`/opt/scRNA_workflow` inside the container. Workflow arguments must refer to paths accessible inside
the container. These checks verify startup and executable availability; they do not run an analysis.

## 5. Installation troubleshooting

| Symptom | What to check |
|---|---|
| `apptainer: command not found` | Install or load the host Apptainer runtime before running the workflow |
| Missing source file during build | Run from the repository root and check the paths in the definition's `%files` section |
| Dependency download failure | Check network access and the first failing download in the build log |
| `No space left on device` | Check the output, temporary, and cache locations, including storage quotas |
| `failed to mount proc` or `operation not permitted` | Check host container permissions and fakeroot support with the system administrator; this error alone does not show that the image is damaged |
| Input file missing inside the container | Check the bind mount and the container-side path |
| Container startup check fails | Resolve the error from `apptainer exec scRNA_seq.sif echo "Container started successfully"` before testing the workflow; analysis parameters cannot fix a container startup failure |


</div>
