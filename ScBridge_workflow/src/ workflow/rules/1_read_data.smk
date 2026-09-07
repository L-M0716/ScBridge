# rules/1_read_data.smk

# -------------Read the sample table to get SAMPLES / RUNS / SAMPLE_RUNS -------------
def _parse_samples_table(path):
    sample_runs = {}
    try:
        with open(path) as f:
            header = next(f).rstrip("\n").split("\t")
    except FileNotFoundError:
        raise ValueError(f"Samples table not found: {path}")
    except StopIteration:
        raise ValueError(f"Samples table is empty: {path}")

    idx = {name: i for i, name in enumerate(header)}
    missing = [c for c in ("sample_id", "run_ids") if c not in idx]
    if missing:
        raise ValueError(f"Samples table missing required columns: {', '.join(missing)}")

    with open(path) as f:
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            sid = fields[idx["sample_id"]]
            runs = [r.strip() for r in fields[idx["run_ids"]].split(",") if r.strip()]
            sample_runs[sid] = runs
    return sample_runs

SAMPLE_RUNS = _parse_samples_table(config["input"]["samples_table"])
SAMPLES = list(SAMPLE_RUNS.keys())
RUNS = sorted({r for runs in SAMPLE_RUNS.values() for r in runs})

def get_sample_runs():
    """return sample_id -> [run_id] Mapping"""
    return SAMPLE_RUNS

def get_all_runs():
    """Return the list of all run_id"""
    return RUNS

def get_runs_for_sample(sample_id):
    """Return the list of run_id for a certain sample"""
    return SAMPLE_RUNS.get(sample_id, [])

wildcard_constraints:
    sample_id="[^/]+",
    run_id="[^/]+"

ruleorder:
    run_starsolo_count > merge_sample_runs


# ------------- checkpoint -------------
checkpoint validate_samples_table:
    input:
        samples_tsv = config["input"]["samples_table"]
    output:
        flag = config["output"]["results"]["base"] + "/read_data/samples_validated.ok",
        sample_ids = config["output"]["results"]["base"] + "/read_data/sample_ids.txt",
        run_ids = config["output"]["results"]["base"] + "/read_data/run_ids.txt"
    log:
        config["output"]["results"]["base"] + "/read_data/logs/samples_validation.log"
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {log}) $(dirname {output.flag})

        if [ ! -s {input.samples_tsv} ]; then
            echo "[ERROR] The sample table is empty or does not exist: {input.samples_tsv}" >&2
            exit 1
        fi
        required_cols=("sample_id" "batch" "group" "run_ids")
        missing_cols=()
        for col in "${{required_cols[@]}}"; do
            if ! head -n 1 {input.samples_tsv} | grep -q "$col"; then
                missing_cols+=("$col")
            fi
        done

        if [ ${{#missing_cols[@]}} -ne 0 ]; then
            echo "[ERROR] The sample table is missing required columns: ${{missing_cols[@]}}" >&2
            exit 2
        fi

        tail -n +2 {input.samples_tsv} | cut -f1 > {output.sample_ids}

        tail -n +2 {input.samples_tsv} | cut -f4 | tr ',' '\n' | sort -u > {output.run_ids}

        sample_count=$(wc -l < {output.sample_ids})
        run_count=$(wc -l < {output.run_ids})

        echo "pass" >> {log}
        echo "Contains $sample_count samples，$run_count 个runs" >> {log}

        touch {output.flag}
        """

rule run_fastqc:
    input:
        samples_valid = config["output"]["results"]["base"] + "/read_data/samples_validated.ok",
        r1 = lambda wildcards: config['input']['fastq_dir'] + "/" + wildcards.run_id + "_1.fastq.gz",
        r2 = lambda wildcards: config['input']['fastq_dir'] + "/" + wildcards.run_id + "_2.fastq.gz"
    output:
        html_r1 = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/{run_id}_1_fastqc.html",
        zip_r1 = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/{run_id}_1_fastqc.zip",
        html_r2 = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/{run_id}_2_fastqc.html",
        zip_r2 = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/{run_id}_2_fastqc.zip"
    params:
        outdir = config["output"]["results"]["base"] + "/read_data/runs/{run_id}"
    log:
        config["output"]["results"]["base"] + "/read_data/runs/{run_id}/logs/fastqc.log"
    threads: config["resources"]["default_threads"]
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}/logs
        echo "Check input file..." >> {log}
        for f in {input.r1} {input.r2}; do
            if [ ! -s "$f" ]; then
                echo "[ERROR] The file does not exist or is empty: $f" >> {log}
                exit 1
            fi
        done

        echo "Start running FastQC..." >> {log}
        fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} >> {log} 2>&1

        for f in {output.html_r1} {output.html_r2}; do
            if [ ! -f "$f" ]; then
                echo "[ERROR] FastQC output file missing: $f" >> {log}
                exit 3
            fi
        done

        echo "FastQC completed" >> {log}
        """


# -------------STARsolo -------------
rule run_starsolo_count:
    input:
        r1 = lambda wildcards: config['input']['fastq_dir'] + "/" + wildcards.run_id + "_1.fastq.gz",
        r2 = lambda wildcards: config['input']['fastq_dir'] + "/" + wildcards.run_id + "_2.fastq.gz",
        samples_valid = config["output"]["results"]["base"] + "/read_data/samples_validated.ok",
        fastqc_r1 = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/{run_id}_1_fastqc.html",
        fastqc_r2 = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/{run_id}_2_fastqc.html",
        genome_file = config["reference"]["paths"]["star_index"] + "/Genome",
        sa_file = config["reference"]["paths"]["star_index"] + "/SA",
        gtf = config["reference"]["paths"]["gtf"],
        cb_whitelist = (config["reference"]["paths"]["whitelists"]["10xv3"]
                        if config["platform"] == "10x_v3"
                        else config["reference"]["paths"]["whitelists"]["bd_v1"])

    output:
        matrix_dir = directory(config["output"]["results"]["base"] + "/read_data/runs/{run_id}/filtered_feature_bc_matrix"),
        matrix_mtx = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/filtered_feature_bc_matrix/matrix.mtx.gz",
        barcodes_tsv = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features_tsv = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/filtered_feature_bc_matrix/features.tsv.gz",
        summary = config["output"]["results"]["base"] + "/read_data/runs/{run_id}/metrics_summary.csv"
    params:
        run_id = "{run_id}",
        outdir = config["output"]["results"]["base"] + "/read_data/runs/{run_id}",
        sjdb_overhang = config["alignment"]["read_length"]-1,
        soloBarcodeReadLength=config["alignment"]["soloBarcodeReadLength"],
        soloUMIlen = config["alignment"]["soloUMIlen"],
        soloCBlen = config["alignment"]["soloCBlen"],
        platform = config["platform"],
        bd_soloCBposition = config["alignment"].get("bd_v1", {}).get("soloCBposition", "0_0_0_8 0_21_0_29 0_43_0_51"),
        bd_soloUMIposition = config["alignment"].get("bd_v1", {}).get("soloUMIposition", "0_52_0_59"),
        soloStrand = config["alignment"]["bd_v1"]["soloStrand"]
    log:
        config["output"]["results"]["base"] + "/read_data/runs/{run_id}/logs/starsolo.log"
    threads: config["resources"]["alignment"]["star"]["threads"]
    resources:
        mem_gb = config["resources"]["alignment"]["star"]["mem_gb"]
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}/logs
        echo "Start running STARsolo for {params.run_id}..." >> {log}

        TEMP_DIR=$(mktemp -d -p {params.outdir})
        echo "Temporary directory: $TEMP_DIR" >> {log}
        STAR_TMP_DIR="$TEMP_DIR/star_tmp"
        if [ "{params.platform}" = "10x_v3" ]; then
            STAR \
                --runThreadN {threads} \
                --genomeDir $(dirname {input.genome_file}) \
                --readFilesIn {input.r2} {input.r1} \
                --readFilesCommand zcat \
                --soloType CB_UMI_Simple \
                --soloCBwhitelist {input.cb_whitelist} \
                --soloUMIlen {params.soloUMIlen} \
                --soloCBlen {params.soloCBlen} \
                --soloBarcodeReadLength {params.soloBarcodeReadLength} \
                --soloFeatures GeneFull \
                --soloCellFilter EmptyDrops_CR \
                --soloOutFileNames solo_ \
                --outFileNamePrefix $TEMP_DIR/starsolo_ \
                --outTmpDir $STAR_TMP_DIR \
                --sjdbGTFfile {input.gtf} \
                --sjdbOverhang {params.sjdb_overhang} \
                --soloCellReadStats Standard \
                --soloStrand Forward \
                --soloUMIdedup 1MM_CR \
                --soloCBmatchWLtype 1MM_multi \
                --soloUMIfiltering MultiGeneUMI_CR \
                --soloMultiMappers Uniform \
                --outSAMtype None \
                --limitBAMsortRAM 50000000000 \
                >> {log} 2>&1
        elif [ "{params.platform}" = "BD_v1" ]; then
            STAR \
                --runThreadN {threads} \
                --genomeDir $(dirname {input.genome_file}) \
                --readFilesIn {input.r2} {input.r1} \
                --readFilesCommand zcat \
                --soloType CB_UMI_Complex \
                --soloCBwhitelist {input.cb_whitelist} \
                --soloCBposition {params.bd_soloCBposition} \
                --soloUMIposition {params.bd_soloUMIposition} \
                --soloBarcodeReadLength 0 \
                --soloFeatures GeneFull \
                --soloCellFilter EmptyDrops_CR \
                --soloOutFileNames solo_ \
                --outFileNamePrefix $TEMP_DIR/starsolo_ \
                --outTmpDir $STAR_TMP_DIR \
                --sjdbGTFfile {input.gtf} \
                --sjdbOverhang {params.sjdb_overhang} \
                --soloCellReadStats Standard \
                --soloStrand {params.soloStrand} \
                --soloUMIdedup 1MM_All \
                --soloCBmatchWLtype 1MM \
                --soloUMIfiltering MultiGeneUMI \
                --soloMultiMappers Uniform \
                --outSAMtype None \
                >> {log} 2>&1
        else
            echo "[ERROR] Unsupported platform type: {params.platform}please use 10x_v3, BD_v1" >> {log}
            exit 1
        fi

        if [ "{params.platform}" = "BD_v1" ]; then
            SOLO_DIR="$TEMP_DIR/starsolo_solo_GeneFull"
        else
            SOLO_DIR="$TEMP_DIR/starsolo_solo_GeneFull"
        fi

        echo "Check STARsolo output..." >> {log}
        if [ ! -d "$SOLO_DIR/filtered" ]; then
            echo "[ERROR] STARsolo did not generate the filtered directory" >> {log}
            ls -la $SOLO_DIR/ >> {log} 2>&1
            exit 1
        fi

        mkdir -p {output.matrix_dir}

        echo "Compress and move the output file..." >> {log}
        gzip -c $SOLO_DIR/filtered/matrix.mtx > {output.matrix_mtx}
        gzip -c $SOLO_DIR/filtered/barcodes.tsv > {output.barcodes_tsv}
        gzip -c $SOLO_DIR/filtered/features.tsv > {output.features_tsv}

        # Generate metrics_summary.csv
        echo "Generate Indicator Summary..." >> {log}
        echo "metric,value" > {output.summary}

        if [ -f "$TEMP_DIR/starsolo_Log.final.out" ]; then
            total_reads=$(grep "Number of input reads" $TEMP_DIR/starsolo_Log.final.out | awk '{{print $NF}}')
            echo "Number of Reads,$total_reads" >> {output.summary}
            if [ -f "$SOLO_DIR/Summary.csv" ]; then
                n_cells=$(tail -n 1 $SOLO_DIR/Summary.csv | cut -d',' -f1)
                echo "Estimated Number of Cells,$n_cells" >> {output.summary}
            fi
        fi

        rm -rf $TEMP_DIR

        echo "STARsolo completed" >> {log}
        """


rule merge_sample_runs:
    input:
        run_matrices = lambda wildcards: [
            config["output"]["results"]["base"] + f"/read_data/runs/{run_id}/filtered_feature_bc_matrix/matrix.mtx.gz"
            for run_id in get_runs_for_sample(wildcards.sample_id)
        ],
        run_barcodes = lambda wildcards: [
            config["output"]["results"]["base"] + f"/read_data/runs/{run_id}/filtered_feature_bc_matrix/barcodes.tsv.gz" 
            for run_id in get_runs_for_sample(wildcards.sample_id)
        ],
        run_features = lambda wildcards: [
            config["output"]["results"]["base"] + f"/read_data/runs/{run_id}/filtered_feature_bc_matrix/features.tsv.gz"
            for run_id in get_runs_for_sample(wildcards.sample_id)
        ]
    output:
        matrix_dir = directory(config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix"),
        matrix_mtx = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix/matrix.mtx.gz",
        barcodes_tsv = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features_tsv = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix/features.tsv.gz"
    params:
        sample_id = "{sample_id}"
    conda:
        config["conda"]["sctest"]
    log:
        config["output"]["results"]["base"] + "/read_data/{sample_id}/logs/merge.log"
    shell:
        """
        set -euo pipefail

        mkdir -p {output.matrix_dir}
        mkdir -p $(dirname {log})

        TEMP_DIR=$(mktemp -d)
        R_SCRIPT="$TEMP_DIR/merge_matrices.R"
        
        cat > "$R_SCRIPT" << 'EOF'
library(Matrix)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]
matrix_files <- args[2:length(args)]

matrices <- list()
for (i in seq_along(matrix_files)) {{
    dir_path <- dirname(matrix_files[i])
    mat <- Read10X(dir_path)
    colnames(mat) <- paste0(colnames(mat), "-", i)
    matrices[[i]] <- mat
}}

merged_mat <- do.call(cbind, matrices)

writeMM(merged_mat, file.path(output_dir, "matrix.mtx"))
write.table(colnames(merged_mat), file.path(output_dir, "barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(rownames(merged_mat), rownames(merged_mat), "Gene Expression"),
            file.path(output_dir, "features.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

system(paste("gzip -f", file.path(output_dir, "matrix.mtx")))
system(paste("gzip -f", file.path(output_dir, "barcodes.tsv")))
system(paste("gzip -f", file.path(output_dir, "features.tsv")))

cat("Merge completed\\n")
cat("Total cell count:", ncol(merged_mat), "\\n")
cat("Total number of genes:", nrow(merged_mat), "\\n")
EOF

        Rscript "$R_SCRIPT" {output.matrix_dir} {input.run_matrices} >> {log} 2>&1

        rm -rf "$TEMP_DIR"

        echo "Matrix merge completed" >> {log}
        """
rule sample_read_complete:
    input:
        matrix_dir = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix",
        matrix_mtx = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix/matrix.mtx.gz",
        barcodes_tsv = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features_tsv = config["output"]["results"]["base"] + "/read_data/{sample_id}/filtered_feature_bc_matrix/features.tsv.gz"
    output:
        flag = config["output"]["results"]["base"] + "/read_data/samples/{sample_id}/read_complete.ok"
    shell:
        "touch {output.flag}"

rule read_data_all:
    input:
        checkpoint_output = lambda wc: checkpoints.validate_samples_table.get().output.flag,
        flags = expand(
            config["output"]["results"]["base"] + "/read_data/samples/{sample_id}/read_complete.ok",
            sample_id=SAMPLES
        )
    output:
        flag = config["output"]["results"]["base"] + "/read_data/all_read_complete.ok"
    shell:
        "touch {output.flag}"
