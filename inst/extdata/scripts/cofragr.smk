# cofragr pipeline

localrules: all, cofrag_cm_merge, cofrag_compartment, compartment_bigwig

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CORE_FACTOR = float(config.get("CORE_FACTOR", 1))
COFRAG_CORE = config.get("COFRAG_CORE", 4)
BOOTSTRAP = config.get("BOOTSTRAP", 50)
STANDARD_COMP = config.get("STANDARD_COMP", "gene_density.hg19")

# Default base directory for data files. Default: ./data
DATA_DIR = config.get("DATA_DIR", os.path.abspath("data"))

# programmatically locate the script 
def find_main_script(script_name):
    import subprocess

    return subprocess.run(
        [
            "Rscript",
            "-e",
            f'cat(paste0(system.file("extdata/scripts/{script_name}", package = "cofragr")))',
        ],
        capture_output=True,
        shell=False,
        text=True,
        check=True,
    ).stdout


rule cofrag_cm:
    input:
        frag="frag/{sid}.frag.bed.gz",
        frag_idx="frag/{sid}.frag.bed.gz.tbi",
    output:
        cm=temp("temp/{sid}.chr{chrom}.cofrag_cm.bed.gz"),
    log: "log/{sid}.chr{chrom}.cofrag_cm.log"
    params:
        label=lambda wildcards: f"cofrag_cm.{wildcards.sid}.chr{wildcards.chrom}",
        bootstrap=BOOTSTRAP,
        main_script=lambda wildcards: find_main_script("cofragr.R"),
    threads: lambda wildcards, attempt: int(COFRAG_CORE * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)
        Rscript {params.main_script} \
        -i {input.frag} \
        -o "$tmpdir" \
        -s {wildcards.sid} \
        -n {threads} \
        --bootstrap {params.bootstrap} \
        --subsample 10000 \
        --seed 1228 \
        --chroms {wildcards.chrom} \
        --min-mapq 30 \
        --min-fraglen 50 --max-fraglen 350 \
        2>&1 | tee {log}

        output_name={wildcards.sid}.cofrag_cm.bed.gz
        mv "$tmpdir"/"$output_name" {output.cm}.tmp
        mv {output.cm}.tmp {output.cm}
        """

rule cofrag_cm_merge:
    input: expand("temp/{{sid}}.chr{chrom}.cofrag_cm.bed.gz", chrom=range(1, 23))
    output: 
        cm="result/{sid}.cofrag_cm.bed.gz",
        cm_index="result/{sid}.cofrag_cm.bed.gz.tbi"
    shell:
        """
        tmpdir=$(mktemp -d)
        output_file="$tmpdir"/output.bed

        idx=0
        for input_file in {input}
        do
            idx=$((idx + 1))
            if [ $idx == 1 ]; then
                zcat "$input_file" > "$output_file"
            else
                zcat "$input_file" | awk 'substr($0,1,1)!="#"' >> "$output_file"
            fi
        done

        bgzip "$output_file"
        tabix -p bed "$output_file".gz

        mv "$output_file".gz {output.cm}.tmp
        mv {output.cm}.tmp {output.cm}
        mv "$output_file".gz.tbi {output.cm_index}        
        """


rule cofrag_compartment:
    input: 
        cm="result/{sid}.cofrag_cm.bed.gz",
        cm_index="result/{sid}.cofrag_cm.bed.gz.tbi"
    output: 
        comp="result/{sid}.compartment.bedGraph.gz",
        comp_index="result/{sid}.compartment.bedGraph.gz.tbi",
    log: "log/{sid}.compartment.log"
    params:
        label=lambda wildcards: f"cofrag_compartment.{wildcards.sid}",
        standard_comp=STANDARD_COMP
    threads: 1
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300
    shell:
        """
        tmpdir=$(mktemp -d)
        cofragr_comp.R \
        -i {input.cm} \
        -o "$tmpdir" \
        -s {wildcards.sid} \
        --res 500000 \
        --standard-compartment {params.standard_comp} \
        2>&1 | tee {log}

        output_name={wildcards.sid}.compartment.bedGraph.gz
        mv "$tmpdir"/"$output_name" {output.comp}.tmp
        mv {output.comp}.tmp {output.comp}
        mv "$tmpdir"/"$output_name".tbi {output.comp_index}
        """


rule compartment_bigwig:
    input: 
        comp="result/{sid}.compartment.bedGraph.gz",
        chrom_sizes="human_g1k_v37.chrom.sizes"
    output:
        "result/{sid}.compartment.bw",
    params:
        label=lambda wildcards: f"compartment_bigwig.{wildcards.sid}",
    threads: 1
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300
    shell:
        """
        tmpdir=$(mktemp -d)

        cat {input.chrom_sizes} |
        bioawk -t '{{$3=$2;$2=0;print}}'|
        bedtools intersect \
        -a {input.comp} -b - |
        sort -k1,1 -k2,2n | bioawk -t '$4!="."' > "$tmpdir"/comp.bedGraph

        bedGraphToBigWig "$tmpdir"/comp.bedGraph {input.chrom_sizes} "$tmpdir"/comp.bw

        mv "$tmpdir"/comp.bw {output}
        """

