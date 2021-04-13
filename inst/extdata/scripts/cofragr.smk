# cofragr pipeline

localrules: all, merge

# >>> Configuration >>>
FULL_CORES = config.get("FULL_CORES", 16)
PART_CORES = config.get("PART_CORES", 4)
MEM_PER_CORE = config.get("MEM_PER_CORE", 1800)
WALL_TIME_MAX = config.get("WALL_TIME_MAX", 2880)
CORE_FACTOR = float(config.get("CORE_FACTOR", 1))
COFRAG_CORE = config.get("COFRAG_CORE", 4)
BOOTSTRAP = config.get("BOOTSTRAP", 50)

# Default base directory for data files. Default: ./data
DATA_DIR = config.get("DATA_DIR", os.path.abspath("data"))


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
    threads: lambda wildcards, attempt: int(COFRAG_CORE * (0.5 + 0.5 * attempt))
    resources:
        mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
        time=WALL_TIME_MAX,
        time_min=300,
        attempt=lambda wildcards, threads, attempt: attempt
    shell:
        """
        tmpdir=$(mktemp -d)
        cofragr.R \
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

# def determine_ifs_cores(chrom, factor=1):
#     if chrom == "X":
#         cores = 6
#     elif chrom == "Y":
#         cores = 3
#     else:
#         chrom = int(chrom)
#         if chrom >= 1 and chrom <= 7:
#             cores = 6
#         elif chrom >= 8 and chrom <= 15:
#             cores = 4
#         elif chrom >= 16 and chrom <= 22:
#             cores = 3
#         else:
#             raise ValueError

#     return round(factor * cores)


# rule ifs:
#     input: 
#         frag=expand("frag/{sid}.hg19.frag.bed.gz", sid=["res_merged", "non_merged"]),
#         gc="data/human_g1k_v37.gc20bp.bed.gz",
#         mappability="data/wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.bedGraph.gz"
#     output: 
#         ifs=temp("results/all.chr{chrom}.ifs.raw.bedGraph.gz"),
#         ifs_idx=temp("results/all.chr{chrom}.ifs.raw.bedGraph.gz.tbi")
#     log: "results/all.chr{chrom}.ifs.raw.log"
#     params:
#         label=lambda wildcards: f"ifs.all.chr{wildcards.chrom}",
#         input_frag=lambda wildcards, input: ":".join(input.frag)
#     threads: lambda wildcards, attempt: round(determine_ifs_cores(wildcards.chrom, CORE_FACTOR) * (1 + 0.5 * attempt))
#     resources:
#         mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
#         time=WALL_TIME_MAX,
#         time_min=300,
#         attempt=lambda wildcards, threads, attempt: attempt
#     shell:
#         """
#         tmpdir=$(mktemp -d)
#         Rscript scripts/cragr.R ifs \
#         --input {params.input_frag} --prefix "$tmpdir"/all.chr{wildcards.chrom} \
#         --gc {input.gc} --mappability {input.mappability} \
#         --chrom {wildcards.chrom} \
#         --min-fraglen 50 --max-fraglen 500 2>&1 | tee {log}

#         output_name=all.chr{wildcards.chrom}.ifs.raw.bedGraph.gz
#         mv "$tmpdir"/"$output_name" {output.ifs}.tmp
#         mv "$tmpdir"/"$output_name".tbi {output.ifs_idx}.tmp
#         mv {output.ifs}.tmp {output.ifs}
#         mv {output.ifs_idx}.tmp {output.ifs_idx}
#         """

# rule hotspot:
#     input: 
#         ifs="results/all.chr{chrom}.ifs.raw.bedGraph.gz",
#         ifs_idx="results/all.chr{chrom}.ifs.raw.bedGraph.gz.tbi",
#         gc="data/human_g1k_v37.gc20bp.bed.gz",
#         mappability="data/wgEncodeDukeMapabilityUniqueness35bp.hs37-1kg.20bp.bedGraph.gz"
#     output: 
#         ifs=temp("results/all.chr{chrom}.ifs.bedGraph.gz"),
#         ifs_idx=temp("results/all.chr{chrom}.ifs.bedGraph.gz.tbi"),
#         hotspot=temp("results/all.chr{chrom}.hotspot.bed.gz"),
#         hotspot_idx=temp("results/all.chr{chrom}.hotspot.bed.gz.tbi"),
#     log: "results/all.chr{chrom}.hotspot.log"
#     params:
#         label=lambda wildcards: f"hotspot.all.chr{wildcards.chrom}",
#     threads: lambda wildcards, attempt: round(determine_ifs_cores(wildcards.chrom, CORE_FACTOR) * (1 + 0.5 * attempt) * 0.667)
#     resources:
#         mem_mb=lambda wildcards, threads: threads * MEM_PER_CORE,
#         time=WALL_TIME_MAX,
#         time_min=300,
#         attempt=lambda wildcards, threads, attempt: attempt
#     shell:
#         """
#         tmpdir=$(mktemp -d)
#         Rscript scripts/cragr.R hotspot \
#         --input {input.ifs} --prefix "$tmpdir"/all.chr{wildcards.chrom} \
#         --gc {input.gc} --mappability {input.mappability} \
#         --chrom {wildcards.chrom} \
#         --min-fraglen 50 --max-fraglen 500 2>&1 | tee {log}

#         output_name=all.chr{wildcards.chrom}.ifs.bedGraph.gz
#         mv "$tmpdir"/"$output_name" {output.ifs}.tmp
#         mv "$tmpdir"/"$output_name".tbi {output.ifs_idx}.tmp
#         mv {output.ifs}.tmp {output.ifs}
#         mv {output.ifs_idx}.tmp {output.ifs_idx}

#         output_name=all.chr{wildcards.chrom}.hotspot.bed.gz
#         mv "$tmpdir"/"$output_name" {output.hotspot}.tmp
#         mv "$tmpdir"/"$output_name".tbi {output.hotspot_idx}.tmp
#         mv {output.hotspot}.tmp {output.hotspot}
#         mv {output.hotspot_idx}.tmp {output.hotspot_idx}
#         """


# rule merge:
#     input: 
#         ifs=expand("results/all.chr{chrom}.ifs.bedGraph.gz", chrom=list(range(1, 23))),
#         ifs_idx=expand("results/all.chr{chrom}.ifs.bedGraph.gz.tbi", chrom=list(range(1, 23))),
#         hotspot=expand("results/all.chr{chrom}.hotspot.bed.gz", chrom=list(range(1, 23))),
#         hotspot_idx=expand("results/all.chr{chrom}.hotspot.bed.gz.tbi", chrom=list(range(1, 23))),
#     output:
#         ifs="results/all.ifs.bedGraph.gz",
#         ifs_idx="results/all.ifs.bedGraph.gz.tbi",
#         hotspot="results/all.hotspot.bed.gz",
#         hotspot_idx="results/all.hotspot.bed.gz.tbi",
#     shell:
#         """
#         tmpdir=$(mktemp -d)

#         set +o pipefail
#         zcat {input.ifs[0]} | head -n 1 | bgzip > $tmpdir/ifs.bedGraph.gz
#         for input_ifs in {input.ifs}
#         do
#             echo Writing $input_ifs ...
#             zcat $input_ifs | tail -n +2 | bgzip >> $tmpdir/ifs.bedGraph.gz
#         done
#         set -o pipefail
#         tabix -p bed $tmpdir/ifs.bedGraph.gz

#         set +o pipefail
#         zcat {input.hotspot[0]} | head -n 1 | bgzip > $tmpdir/hotspot.bed.gz
#         for input_hotspot in {input.hotspot}
#         do
#             echo Writing $input_hotspot ...
#             zcat $input_hotspot | tail -n +2 | bgzip >> $tmpdir/hotspot.bed.gz
#         done
#         set -o pipefail
#         tabix -p bed $tmpdir/hotspot.bed.gz

#         mv "$tmpdir"/ifs.bedGraph.gz {output.ifs}.tmp
#         mv "$tmpdir"/ifs.bedGraph.gz.tbi {output.ifs_idx}.tmp
#         mv {output.ifs}.tmp {output.ifs}
#         mv {output.ifs_idx}.tmp {output.ifs_idx}

#         mv "$tmpdir"/hotspot.bed.gz {output.hotspot}.tmp
#         mv "$tmpdir"/hotspot.bed.gz.tbi {output.hotspot_idx}.tmp
#         mv {output.hotspot}.tmp {output.hotspot}
#         mv {output.hotspot_idx}.tmp {output.hotspot_idx}
#         """

# sample_ids, = glob_wildcards("frag/{sample,Pilot2.+}.hg19.frag.bed.gz")
# # Head & neck cancer samples
# responder_ids = [12, 25, 36, 39, 68, 78, 80, 82]
# non_responder_ids = [1, 3, 6, 9, 14, 16, 22, 28, 30 ,34, 42, 45, 63, 65, 71, 75, 86, 89]
# sample_ids = responder_ids + non_responder_ids

# rule all:
#     input: expand("results/Pilot2-{sample}.ifs.bedGraph.gz", sample = sample_ids)