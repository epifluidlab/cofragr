# cofragr pipeline

localrules: all, cofrag_cm_merge, cofrag_compartment, compartment_bigwig

# >>> Configuration >>>
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)

# cofragr specific configuration
COFRAG_CORE = config.get("COFRAG_CORE", 4)
SUBSAMPLE = config.get("SUBSAMPLE", 10000)
MIN_MAPQ = config.get("MIN_MAPQ", 30)
MIN_FRAGLEN = config.get("MIN_FRAGLEN", 50)
MAX_FRAGLEN = config.get("MAX_FRAGLEN", 350)
BOOTSTRAP = config.get("BOOTSTRAP", 50)
SEED = config.get("SEED", 1228)
BLOCK_SIZE = config.get("BLOCK_SIZE", 10000000)
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
        subsample=SUBSAMPLE,
        seed=SEED,
        min_mapq=MIN_MAPQ,
        min_fraglen=MIN_FRAGLEN,
        max_fraglen=MAX_FRAGLEN,
        block_size=BLOCK_SIZE,
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
        --block-size {params.block_size} \
        --bootstrap {params.bootstrap} \
        --subsample {params.subsample} \
        --seed {params.seed} \
        --chroms {wildcards.chrom} \
        --min-mapq {params.min_mapq} \
        --min-fraglen {params.min_fraglen} \
        --max-fraglen {params.max_fraglen} \
        2>&1 | tee {log}

        ls -la "$tmpdir"
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

