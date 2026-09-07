# rules/0_preprocess.smk

### rule 1: Reference Data Basic Verification
rule check_reference_data:
    input:
        genome_fa = config["reference"]["paths"]["genome_fa"],
        genome_gtf = config["reference"]["paths"]["gtf"]
    output:
        flag = config["output"]["results"]["base"] + "/preprocess/ref_check.ok"
    log:
        config["output"]["results"]["base"] + "/preprocess/logs/reference_validation.log"
    shell:
        """
        set -euo pipefail
        
        mkdir -p $(dirname {log}) $(dirname {output.flag})
        
        missing_files=()
        for f in {input}; do
            if [ ! -f "$f" ]; then
                missing_files+=("$f")
            fi
        done
        
        if [ ${{#missing_files[@]}} -ne 0 ]; then
            echo "[ERROR] Missing reference document: ${{missing_files[@]}}" >&2
            echo "[ERROR] Missing reference document: ${{missing_files[@]}}" >> {log}
            exit 1
        fi
        
        # 1. FASTA file
        if ! head -n 1 {input.genome_fa} | grep -q '^>'; then
            echo "[ERROR] Genomic FASTA format abnormal" >&2
            exit 2
        fi
        
        if [ $(wc -l < {input.genome_fa}) -lt 2 ]; then
            echo "[ERROR] The genome FASTA file is incomplete" >&2
            exit 3
        fi
        
        if ! grep -m 1 -q 'gene_id' {input.genome_gtf}; then
            echo "[ERROR] The GTF file is missing the gene_id field" >&2
            exit 4
        fi
        echo "pass ($(date))" >> {log}
        touch {output.flag}
        """

rule build_star_index:
    input:
        fa = config["reference"]["paths"]["genome_fa"],
        gtf = config["reference"]["paths"]["gtf"],
        ref_check = config["output"]["results"]["base"] + "/preprocess/ref_check.ok"
    output:
        index_flag = config["output"]["results"]["base"] + "/preprocess/star_index.ok",
        genome_file = config["reference"]["paths"]["star_index"] + "/Genome",
        sa_file = config["reference"]["paths"]["star_index"] + "/SA",
        params_file = config["reference"]["paths"]["star_index"] + "/genomeParameters.txt"
    params:
        idx_dir = config["reference"]["paths"]["star_index"],
        num_threads = config["resources"]["alignment"]["star"]["index_threads"],
        sjdb_overhang = config["alignment"]["read_length"] - 1
    log:
        config["output"]["results"]["base"] + "/preprocess/logs/star_index.log"
    shell:
        """
        set -euo pipefail
        
        mkdir -p $(dirname {log}) {params.idx_dir}
        
        if [ -f "{output.genome_file}" ] && [ -f "{output.sa_file}" ]; then
            echo "pass build" > {log}
        else
            echo "=== STAR index construction started ($(date)) ===" > {log}
            echo "Parameter:" >> {log}
            echo "  Number of threads: {params.num_threads}" >> {log}
            echo "  sjdbOverhang: {params.sjdb_overhang}" >> {log}
            
            STAR \
                --runMode genomeGenerate \
                --genomeDir {params.idx_dir} \
                --genomeFastaFiles {input.fa} \
                --sjdbGTFfile {input.gtf} \
                --runThreadN {params.num_threads} \
                --sjdbOverhang {params.sjdb_overhang} >> {log} 2>&1
            for f in {output.genome_file} {output.sa_file} {output.params_file}; do
                if [ ! -f "$f" ]; then
                    echo "error: $f" >> {log}
                    exit 1
                fi
            done
            echo "Index construction completed" >> {log}
        fi
        
        touch {output.index_flag}
        """
rule preprocess_all:
    input:
        config["output"]["results"]["base"] + "/preprocess/ref_check.ok",
        config["output"]["results"]["base"] + "/preprocess/star_index.ok"
