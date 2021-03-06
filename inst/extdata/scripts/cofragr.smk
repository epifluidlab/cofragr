# cofragr pipeline
from snakemake.utils import min_version
min_version("6.0")

# localrules: all, cofrag_cm_merge, cofrag_compartment, compartment_bigwig

# >>> Configuration >>>
MEM_PER_CORE = config.get("MEM_PER_CORE", 2000)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
WALL_TIME_MIN = config.get("WALL_TIME_MIN", 300)

# cofragr specific configuration
# COFRAG_MEM = int(config.get("COFRAG_MEM"), 6000)
COFRAG_CORE = int(config.get("COFRAG_CORE", 4))

SUBSAMPLE = int(config.get("SUBSAMPLE", 10000))
MIN_MAPQ = int(config.get("MIN_MAPQ", 30))
MIN_FRAGLEN = int(config.get("MIN_FRAGLEN", 50))
MAX_FRAGLEN = int(config.get("MAX_FRAGLEN", 350))
BOOTSTRAP = int(config.get("BOOTSTRAP", 50))
SEED = int(config.get("SEED", 1228))
BLOCK_SIZE = int(config.get("BLOCK_SIZE", 10000000))
BIN_SIZE = int(config.get("BIN_SIZE", 500000))

CHROM_LIST = config.get("CHROM_LIST", ",".join([str(v) for v in range(1, 23)] + ["X"])).split(",")
# STANDARD_COMP = config.get("STANDARD_COMP", "gene_density.hg19")
STANDARD_COMP = config.get("STANDARD_COMP", "wbc")

# Default base directory for data files. Default: ./data
DATA_DIR = config.get("DATA_DIR", os.path.abspath("data"))

SCRIPT_PATH = config.get("SCRIPT_PATH", ".") 
DISABLE_PARALLEL = bool(config.get("DISABLE_PARALLEL", 0))
R = config.get("R", "Rscript")


def mem_for_attempt(mem_mb, attempt):
    return int(mem_mb * (0.5 + 0.5 * attempt))


def threads_for_mem(mem_mb, mem_per_core=MEM_PER_CORE):
    return int(math.ceil(mem_mb / mem_per_core) + 1)


rule cofrag_cm:
    input:
        frag="frag/{sid}.{genome}.frag.bed.gz",
        frag_idx="frag/{sid}.{genome}.frag.bed.gz.tbi",
    output:
        cm=temp("temp/{sid}.{genome,(hg19|hg38|GRCh37|GRCh38)}.chr{chrom}.cofrag_cm.bed.gz"),
    log: "log/{sid}.{genome}.chr{chrom}.cofrag_cm.log"
    params:
        slurm_job_label=lambda wildcards: f"cofrag_cm.{wildcards.sid}.{wildcards.genome}.chr{wildcards.chrom}",
        bootstrap=BOOTSTRAP,
        subsample=SUBSAMPLE,
        seed=SEED,
        min_mapq=MIN_MAPQ,
        min_fraglen=MIN_FRAGLEN,
        max_fraglen=MAX_FRAGLEN,
        block_size=BLOCK_SIZE,
        bin_size=BIN_SIZE,
        disable_parallel=lambda wildcards: "--disable-parallel" if DISABLE_PARALLEL else "",
    threads: lambda wildcards, attempt: int(COFRAG_CORE * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        {R} -e 'sessionInfo()'

        {R} {SCRIPT_PATH}/cofragr.R \
        -i {input.frag} \
        -o "$tmpdir" \
        -s {wildcards.sid} \
        -g {wildcards.genome} \
        -n {threads} \
        --res {params.bin_size} \
        --block-size {params.block_size} \
        --bootstrap {params.bootstrap} \
        --subsample {params.subsample} \
        --seed {params.seed} \
        --chroms {wildcards.chrom} \
        --min-mapq {params.min_mapq} \
        --min-fraglen {params.min_fraglen} \
        --max-fraglen {params.max_fraglen} \
        {params.disable_parallel} \
        2>&1 | tee {log}

        output_name={wildcards.sid}.cofrag_cm.bed.gz
        mv $tmpdir/$output_name {output.cm}.tmp
        mv {output.cm}.tmp {output.cm}

        rm -rf $tmpdir
        """

rule cofrag_cm_merge:
    input: expand("temp/{{sid}}.{{genome}}.chr{chrom}.cofrag_cm.bed.gz", chrom=CHROM_LIST)
    output:
        cm="cofrag/{sid}.{genome,(hg19|hg38|GRCh37|GRCh38)}.cofrag_cm.bed.gz",
        cm_index="cofrag/{sid}.{genome}.cofrag_cm.bed.gz.tbi"
    threads: lambda wildcards, attempt: int(2 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    params:
        slurm_job_label=lambda wildcards: f"cofrag_cm_merge.{wildcards.sid},{wildcards.genome}",
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        output_file=$tmpdir/output.bed.gz

        idx=0
        for input_file in {input}
        do
            echo "Processing $input_file" ...
            idx=$((idx + 1))
            if [ $idx == 1 ]; then
                cp "$input_file" "$output_file"
            else
                zcat "$input_file" | awk 'substr($0,1,1)!="#"' | bgzip >> "$output_file"
            fi
        done

        tabix -p bed "$output_file"

        mv "$output_file" {output.cm}.tmp
        mv {output.cm}.tmp {output.cm}
        mv "$output_file".tbi {output.cm_index}

        rm -rf $tmpdir
        """


rule cofrag_compartment:
    input:
        cm="cofrag/{sid}.{genome}.cofrag_cm.bed.gz",
        cm_index="cofrag/{sid}.{genome}.cofrag_cm.bed.gz.tbi"
    output:
        comp="cofrag/{sid}.{genome,(hg19|hg38|GRCh37|GRCh38)}.compartment.{method,(juicer|lieberman|obs_exp)}.bed.gz",
        #comp_correlation="cofrag/{sid}.{genome}.compartment_correlation.{method}.tsv",
    log: "log/{sid}.{genome}.compartment.{method}.log"
    params:
        slurm_job_label=lambda wildcards: f"cofrag_compartment.{wildcards.sid}.{wildcards.genome}.{wildcards.method}",
        standard_comp=STANDARD_COMP,
        bin_size=BIN_SIZE,
    threads: lambda wildcards, attempt: int(8 * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        set +u; if [ -z $LOCAL ] || [ -z $SLURM_CLUSTER_NAME ]; then tmpdir=$(mktemp -d); else tmpdir=$(mktemp -d -p $LOCAL); fi; set -u

        {R} -e 'sessionInfo()'

        {R} {SCRIPT_PATH}/cofragr_comp.R \
        -o $tmpdir/output.bed.gz \
        --res {params.bin_size} \
        --method {wildcards.method} \
        --genome {wildcards.genome} \
        --reference gc \
        --smooth 1,2,3 \
        -n {threads} \
        {input.cm} \
        2>&1 | tee {log}

        mv $tmpdir/output.bed.gz {output.comp}.tmp
        mv {output.comp}.tmp {output.comp}
        """


rule compartment_bigwig:
    input:
        comp="cofrag/{sid}.compartment.bedGraph.gz",
        chrom_sizes="human_g1k_v37.chrom.sizes"
    output:
        "cofrag/{sid}.compartment.bw",
    params:
        slurm_job_label=lambda wildcards: f"compartment_bigwig.{wildcards.sid}",
    threads: 1
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=WALL_TIME_MIN
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

