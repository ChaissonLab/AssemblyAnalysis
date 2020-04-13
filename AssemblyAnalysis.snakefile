import os
import tempfile
import subprocess
import os.path
env = os.environ
import sys
if "TMPDIR" not in env:
    print("ERROR, TMPDIR must be defined")
    sys.exit(1)

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "assembly_analysis.json"



haps=["1", "2"]
ref=config["ref"]

fai = open(ref + ".fai")
allChroms = [l.strip().split()[0] for l in fai ]
chroms = []

asms=list(config["assemblies"].keys())

# Index assemblies by name and haplotype
asmFiles={}
for asm in config["assemblies"]:
    asmFiles[asm] = {}
    asmFiles[asm]["assemb"] = {}
    asmFiles[asm]["assemb"]["1"] = config["assemblies"][asm][0]
    asmFiles[asm]["assemb"]["2"] = config["assemblies"][asm][1]
    asmFiles[asm]["readType"] = config["assemblies"][asm][2]

methods=["lra", "mm2"]    
ops=["ins", "del"]
cats=["missed", "found"]
cats +=["missed_slop", "found_slop"]
checks=["passed","failed"]
flags=["1", "2", "3"] #Current number of filtering flags

localrules:IndexAssembly,FindAvgCov,MakeBed,SplitBed,GetVariantSupport,FlagCentromere,FlagUnsupported,FlagHighCoverage,CollectFlagStats,CollectRunStats,CombineRunStats,CombineFlagStats,CompareToSVSet,CombineVariantSupport,CombineChroms,ReadVCFToBed,ClusterBed,ClusterDiploid

rule all:
    input:
        readVCF=expand("from_reads.{method}.vcf",method=methods),
        readBed=expand("from_reads.{method}.{op}.bed",method=methods,op=ops),
        aligned=expand("{asm}.{hap}/{method}/aln.bam",asm=asms,hap=haps,method=methods),
        asmBed=expand("{asm}.{hap}/{method}/variants.sv.bed",asm=asms,hap=haps,method=methods),
        splitBed=expand("{asm}.{hap}/{method}/variants.sv.{op}.bed",asm=asms,hap=haps,method=methods,op=ops),
        supBed=expand("{asm}.{hap}/{method}/variants.sv.{op}.bed.{meth}.support",asm=asms,hap=haps,method=methods,op=ops,meth=methods),
        supComb=expand("{asm}.{hap}/{method}/variants.sv.{op}.bed.support.combined",asm=asms,hap=haps,method=methods,op=ops),
        clustBed=expand("{asm}.{hap}/{method}/variants.clusters.bed",asm=asms,hap=haps,method=methods),
        clustDip=expand("{asm}_dip/{method}/variants.clustered.bed",asm=asms,method=methods),
        comparison=expand("{asm}_dip/{method}/variants.sv.{op}.{cat}.bed",asm=asms,method=methods,op=ops+["hap"],cat=cats),
        filtered=expand("{asm}.{hap}/{method}/sv.{check}.bed",asm=asms,hap=haps,method=methods,check=checks),
        findCov=expand("coverage/{method}.txt",method=methods),
        statsAsm=expand("{asm}.{hap}/{method}/stats.txt",asm=asms,hap=haps,method=methods),
        statsAsmFlags=expand("{asm}.{hap}/{method}/stats_flag.txt",asm=asms,hap=haps,method=methods),
        statsAll="run_stats_all.txt",
        statsAllFlags="run_stats_all_flags.txt"


rule IndexAssembly:
    input:
        asm=lambda wildcards: asmFiles[wildcards.asm]["assemb"][wildcards.hap]
    output:
        index="indices/{asm}.{hap}.assembly.fasta.fai"
    shell:"""
mkdir -p indices
ln -sfn {input.asm} indices/{wildcards.asm}.{wildcards.hap}.assembly.fasta && samtools faidx indices/{wildcards.asm}.{wildcards.hap}.assembly.fasta && rm -f indices/{wildcards.asm}.{wildcards.hap}.assembly.fasta
"""

rule FindAvgCov:
    input:
        reads=lambda wildcards: config["bam"][wildcards.method]
    output:
        cov="coverage/{method}.txt"
    shell:"""
mkdir -p coverage
samtools coverage -H {input.reads} > {output.cov}
"""

#
#Step 1: Map contigs
#
rule MapFasta:
    input:
        asm=lambda wildcards: asmFiles[wildcards.asm]["assemb"][wildcards.hap]
    output:
        bam="{asm}.{hap}/{method}/aln.bam"
    resources:
        threads=4
    params:
        grid_opts=config["grid_memory"],
        readType=lambda wildcards: asmFiles[wildcards.asm]["readType"],
        ref=config["ref"]
    shell:"""
mkdir -p {wildcards.asm}.{wildcards.hap}/{wildcards.method}

if [ {wildcards.method} == mm2 ] ; then \
   minimap2 -t 4 -ax asm10 --cs {params.ref} {input.asm} | samtools view -h -b - ; \
elif [ {wildcards.method} == lra ] ; then \
   lra align -{params.readType} -t 4 {params.ref} {input.asm} -p s ; fi | samtools sort -@2 -T $TMPDIR/$$  -o {output.bam}

samtools index {output.bam}
"""

#
#Step 2: Call SVs (variants >= 50bp) from mapped assemblies. Note: For minimap2, Heng's paftools pipeline is used.
#
rule MakeBed:
    input:
        bam="{asm}.{hap}/{method}/aln.bam"
    output:
        bed="{asm}.{hap}/{method}/variants.sv.bed"
    params:
        ref=config["ref"],
        sd=SD
    shell:"""

samtools view -h {input.bam} | \
   if [ {wildcards.method} == mm2 ] ; then \
     paftools.js sam2paf - | sort -k6,6 -k8,8n | paftools.js call -f {params.ref} - | \
     grep -v "^#" | awk '{{diff=length($4)-length($5); if (diff <= -50 || diff >= 50) {{ \
                         if (diff <= 0) {{oper="insertion"; diff*= -1}} else {{oper="deletion"}}; \
                         print $1,$2,$2+diff,oper,diff,$4,$5}} }}' OFS="\t"; \
   elif [ {wildcards.method} == lra ] ; then \
     {params.sd}/PrintGaps.py {params.ref} /dev/stdin | tail -n +2 | cut -f1-7 ; fi |\
   awk -F " " '!keep[$1,$2,$3,$4]++' OFS="\t" | bedtools sort > {output.bed}
"""

#
#Step 3: Create bed files with SV calls split by operation: deletion and insertion.
#
rule SplitBed:
    input:
        bed="{asm}.{hap}/{meth}/variants.sv.bed"
    output:
        split=expand("{{asm}}.{{hap}}/{{meth}}/variants.sv.{op}.bed",op=ops)
    shell:"""
beds[0]={output.split[0]} ; beds[1]={output.split[1]}
oper=(insertion deletion) ; for op in {{0..1}} ; do egrep "^#|${{oper[$op]}}" {input.bed} | \
       awk 'function max(x, y) {{return length(x) > length(y) ? x: y}} {{print $1,$2,$3,$4,$5,max($6,$7)}}' OFS="\t" > ${{beds[$op]}}; done
"""

#
#Step 4: Compare SV calls to existing diploid SV set.
#
rule ClusterBed:
    input:
        bed="{asm}.{hap}/{method}/variants.sv.bed"
    output:
        clustBed="{asm}.{hap}/{method}/variants.clusters.bed"
    shell:"""
bedtools cluster -d 500 -i {input.bed} | bedtools groupby -g 8 -c 1,2,3,4,5,6,7,8 -o first,min,max,collapse,collapse,collapse,collapse,count | cut -f 2- | \
awk '{{ if (index($4,"deletion") > 0 && index($4,"insertion") > 0) {{ $4="locus";}} else {{ if (index($4,"deletion") == 0) {{ $4="insertion"; }} else {{ $4="deletion";}} }} print; }}' | tr " " "\\t"  > {output.clustBed}
"""

rule ClusterDiploid:
    input:
        clustBed=expand("{{asm}}.{hap}/{{method}}/variants.clusters.bed",hap=haps)
    output:
        dip="{asm}_dip/{method}/variants.clustered.bed",
    shell:"""
cat {input.clustBed} | bedtools sort | bedtools cluster -d 500 | bedtools groupby -c 9 -g 9 -o count -full > {output.dip}
"""

rule CompareToSVSet:
    input:
        bed="{asm}_dip/{method}/variants.clustered.bed",
        set= lambda wildcards: config["sv-set"][wildcards.op]
    output:
        missed="{asm}_dip/{method}/variants.sv.{op}.missed.bed",
        found="{asm}_dip/{method}/variants.sv.{op}.found.bed",
        missed_slop="{asm}_dip/{method}/variants.sv.{op}.missed_slop.bed",
        found_slop="{asm}_dip/{method}/variants.sv.{op}.found_slop.bed"
    params:
        ref=config["ref"]
    shell:"""

bedtools intersect -v -a {input.set} -b {input.bed} > {output.missed}
bedtools intersect -u -a {input.set} -b {input.bed} > {output.found}

bedtools slop -i {input.bed} -b 1000 -g {params.ref}.fai | bedtools intersect -v -a {input.set} -b stdin > {output.missed_slop}
bedtools slop -i {input.bed} -b 1000 -g {params.ref}.fai | bedtools intersect -u -a {input.set} -b stdin > {output.found_slop}
"""

#
#Step 5: Create set of SVs (variants >= 50bp) called from raw reads.
#
rule CreateChromVCF:
    input:
        reads=lambda wildcards: config["bam"][wildcards.method]
    output:
        vcf="gaps/{method}.{chrom}.gaps.vcf"
    resources:
        threads=16
    params:
        grid_opts=config["grid_large"],
        ref=config["ref"],
        sample=config["sample"],
        sd=SD
    shell:"""
mkdir -p gaps
{params.sd}/SamToVCF.py --ref {params.ref} --sample {params.sample} --sam {input.reads} --minLength 25 --chrom {wildcards.chrom} > {output.vcf}
"""

rule CombineChroms:
    input:
        vcfs=expand("gaps/{{method}}.{chrom}.gaps.vcf", chrom=allChroms)
    output:
        gapVcf="from_reads.{method}.vcf"
    shell:"""
grep "^#" {input.vcfs[0]} > {output.gapVcf}
cat {input.vcfs} | grep -v "^#" >> {output.gapVcf}
"""

rule ReadVCFToBed:
    input:
        vcf="from_reads.{method}.vcf"
    output:
        bed=expand("from_reads.{{method}}.{op}.bed",op=ops)
    params:
        sd=SD
    shell:"""
{params.sd}/VcfToSplitBed.sh {input.vcf} {output.bed[1]} {output.bed[0]}
"""

#
#Step 6: Compare SVs called from the assembly to SVs called from the raw reads to determine if the assembly-SVs are supported..
#
rule GetVariantSupport:
    input:
        bed="{asm}.{hap}/{method}/variants.sv.{op}.bed",
        sup="from_reads.{method}.{op}.bed"
    output:
        sup="{asm}.{hap}/{method}/variants.sv.{op}.bed.{meth}.support"
    params:
        grid_opts=config["grid_medium"],
        ref=config["ref"],
        sd=SD
    shell:"""
bedtools slop -g {params.ref}.fai -b 1000 -i {input.sup} | \
  bedtools intersect -loj -a {input.bed} -b stdin | \
  bedtools groupby -c 11 -o collapse -full | \
  {params.sd}/GetSVSupport.py > {output.sup}
"""

rule CombineVariantSupport:
    input:
        meth=expand("{{asm}}.{{hap}}/{{method}}/variants.sv.{{op}}.bed.{meth}.support", meth=methods),
    output:
        comb="{asm}.{hap}/{method}/variants.sv.{op}.bed.support.combined"
    shell:"""
paste {input.meth[0]} <( cut -f 7 {input.meth[1]} ) > {output.comb}
"""

#
#Step 7: Filter false-positives out of the set of SVs called from the assembly
#
rule FlagCentromere:
    input:
        sup=expand("{{asm}}.{{hap}}/{{meth}}/variants.sv.{op}.bed.support.combined",op=ops)
    output:
        euchr="{asm}.{hap}/{meth}/sv.flag.1.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""

bedtools intersect -c -a <(cat {input.sup}) -b /home/cmb-16/mjc/shared/references/hg38/regions/cytobands/hg38.heterochromatic.bed | awk '{{if ($9 > 0) print $0,1 ; else print $0,0}}' OFS="\t" > {output.euchr}
"""

rule FlagHighCoverage:
    input:
        sup=expand("{{asm}}.{{hap}}/{{meth}}/variants.sv.{op}.bed.support.combined",op=ops),
        reads=lambda wildcards: config["bam"][wildcards.meth],
        cov="coverage/{meth}.txt",
    output:
        hcov="{asm}.{hap}/{meth}/sv.flag.2.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""

avgCov=$(cat {input.cov} | awk '{{ total+=$7 }} END {{ print total/NR }}')

samtools bedcov <(cat {input.sup}) {input.reads} | awk -v avg=$avgCov '{{if ($9/$5 > 2*avg) print $0,1 ; else print $0,0}}' OFS="\t" > {output.hcov}
"""

rule FlagUnsupported:
    input:
        sup=expand("{{asm}}.{{hap}}/{{meth}}/variants.sv.{op}.bed.support.combined",op=ops)
    output:
        unsup="{asm}.{hap}/{meth}/sv.flag.3.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""

cat {input.sup} | awk '{{if ($7 < 3 && $8 < 3) print $0,1 ; else print $0,0}}' OFS="\t" > {output.unsup}
"""

rule CheckFlags:
    input:
        flagged=expand("{{asm}}.{{hap}}/{{meth}}/sv.flag.{flag}.bed",flag=flags)
    output:
        filtered=expand("{{asm}}.{{hap}}/{{meth}}/sv.{check}.bed",check=checks)
    params:
        grid_opts=config["grid_small"]
    shell:"""

touch {output.filtered}

paste {input.flagged[0]} <( cut -f 10 {input.flagged[1]} ) <( cut -f 9 {input.flagged[2]} ) | awk -v asm="{wildcards.asm}" -v meth="{wildcards.meth}" -v hap="{wildcards.hap}" '{{if ($10+$11+$12 == 0) print $1,$2,$3,asm"."hap"/"meth"/"$4,$5,$6 > "{output.filtered[0]}" ; else print $1,$2,$3,asm"."hap"/"meth"/"$4,$5,$6,$10,$11,$12 > "{output.filtered[1]}"; }}' OFS="\t"
"""

#
#Step 8: Collect statistics on this run of the pipeline
#
rule CollectRunStats:
    input:
        asm="indices/{asm}.{hap}.assembly.fasta.fai",
        covs="{asm}.{hap}/{method}/aln.bam"
    output:
        stats="{asm}.{hap}/{method}/stats.txt"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"]
    shell:"""

halfGnm=$(cat {params.ref}.fai | awk '{{size += $2}} END {{print int(size/2)}}')

n50=$(sort -nr -k 2 {input.asm} | awk -v gnm="$halfGnm" '{{contgSum += $2 ; if(contgSum >= gnm) {{print $2}} }} END{{if(contgSum < gnm) {{print "notFound"}} }}' > {output.stats}.tmp ; head -n 1 {output.stats}.tmp) && rm -f {output.stats}.tmp

cov=$(samtools coverage -H {input.covs} | awk '{{count+=$5; total+=$3}} END {{print count/total}}')

echo -e "{wildcards.method}.{wildcards.asm}.{wildcards.hap}\t$n50\t$cov" >> {output.stats}
"""

rule CollectFlagStats:
    input:
        sup=expand("{{asm}}.{{hap}}/{{method}}/sv.{check}.bed",check=checks)
    output:
        flags="{asm}.{hap}/{method}/stats_flag.txt"
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"]
    shell:"""

pass=$(cat {input.sup[0]} | wc -l)

fail=$(cat {input.sup[1]} | wc -l) ; centro=$(cat {input.sup[1]} | awk '{{if ($7==1 && $8+$9==0) print}}' | wc -l) 
high=$(cat {input.sup[1]} | awk '{{if ($8==1 && $7+$9==0) print}}' | wc -l) ; unsup=$(cat {input.sup[1]} | awk '{{if ($9==1 && $7+$8==0) print}}' | wc -l)
centro_high=$(cat {input.sup[1]} | awk '{{if ($7+$8==2 && $9==0) print}}' | wc -l) ; centro_unsup=$(cat {input.sup[1]} | awk '{{if ($7+$9==2 && $8==0) print}}' | wc -l)
high_unsup=$(cat {input.sup[1]} | awk '{{if ($8+$9==2 && $7==0) print}}' | wc -l) ; centro_high_unsup=$(cat {input.sup[1]} | awk '{{if ($7+$8+$9==3) print}}' | wc -l)

total=$(($pass+$fail)) ; total_centro=$(($total-$centro-$centro_high-$centro_unsup-$centro_high_unsup)) ; total_centro_high=$(($total_centro-$high-$high_unsup)) ; full_sup=$(($total_centro_high-$unsup))

echo -e "{wildcards.method}.{wildcards.asm}.{wildcards.hap}\t$total\t$total_centro\t$total_centro_high\t$full_sup\t$centro\t \
         $high\t$unsup\t$centro_high\t$centro_unsup\t$high_unsup\t$centro_high_unsup" >> {output.flags}
"""

rule CombineRunStats:
    input:
        stats=expand("{asm}.{hap}/{method}/stats.txt",asm=asms,hap=haps,method=methods),
        flags=expand("{asm}.{hap}/{method}/stats_flag.txt",asm=asms,hap=haps,method=methods)
    output:
        all="run_stats_all.txt"
    shell:"""

echo -e "Assembly\tN50\tBasesCov\tTotal\tTotal-Centro\tTotal-Centro-HighCov\tSup\n" | cat - <(paste <(cat {input.stats}) <(cat {input.flags} | cut -f 2-5)) > {output.all} 
"""

rule CombineFlagStats:
    input:
        stats=expand("{asm}.{hap}/{method}/stats_flag.txt",asm=asms,hap=haps,method=methods)
    output:
        all="run_stats_all_flags.txt"
    shell:"""

echo -e "Assembly\tTotal\tTotal-Centro\tTotal-Centro-HighCov\tFullSup\tCentro\tHighCov\tUnSup\tCentroHigh\tCentroUnsup\tHighUnsup\tCentroHighUnsup" | cat - {input.stats} > {output.all}
"""
