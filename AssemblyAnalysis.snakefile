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


localrules:IndexAssembly,FindAvgCov,MakeBed,SplitBed,CombineStats,CompareToSVSet,CombineVariantSupport,CombineChroms,ReadVCFToBed,ClusterBed,ClusterDiploid

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
        statsAsm=expand("{asm}.{hap}/stats.txt",asm=asms,hap=haps),
        statsAll="run_stats_all.txt"


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
samtools coverage {input.reads} > {output.cov}
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

if [ {wildcards.method} == "mm2" ] ; then \
   minimap2 -t 4 -ax asm5 --cs {params.ref} {input.asm} | samtools view -h -b - ; \
elif [ {wildcards.method} == "lra" ] ; then \
   /home/cmb-16/mjc/mchaisso/projects/LRA/lra align {params.ref} {input.asm} -{params.readType} -t 4 -p s ; fi | \
   samtools sort -@2 -T $TMPDIR/$$  -o {output.bam}

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
   if [ {wildcards.method} == "mm2" ] ; then \
     paftools.js sam2paf - | sort -k6,6 -k8,8n | paftools.js call -f {params.ref} - | \
     grep -v "^#" | awk '{{diff=length($4)-length($5); if (diff <= -50 || diff >= 50) {{ \
                         if (diff <= 0) {{oper="insertion"; diff*= -1}} else {{oper="deletion"}}; \
                         print $1,$2,$2+diff,oper,diff,$4,$5}} }}' OFS="\t"; \
   elif [ {wildcards.method} == "lra" ] ; then \
     {params.sd}/PrintGaps.py {params.ref} /dev/stdin | tail -n +2 | cut -f1-7 ; fi > {output.bed}
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
       awk 'function max(x, y) {{return length(x) > length(y) ? x: y}} {{print $1,$2,$3,$4,$5,max($6,$7); \
       !keep[$1,$2,$3,$4]++}}' OFS="\t" | bedtools sort > ${{beds[$op]}}; done
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
rule CheckResults:
    input:
        sup=expand("{{asm}}.{{hap}}/{{meth}}/variants.sv.{op}.bed.support.combined",op=ops),
        reads=lambda wildcards: config["bam"][wildcards.meth],
        cov="coverage/{meth}.txt",
    output:
        failed="{asm}.{hap}/{meth}/sv.failed.bed",
        passed="{asm}.{hap}/{meth}/sv.passed.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""

avgCov=$(cat {input.cov} | awk '{{ total += $7 }} END {{ print total/NR }}')

cat {input.sup} | awk -v asm="{wildcards.asm}" -v meth="{wildcards.meth}" -v hap="{wildcards.hap}" '{{if ($7 < 3 && $8 < 3) {{print $1,$2,$3,asm"."hap"/"meth"/"$4,$5,$6,$7,$8 }} }}' OFS="\t" > {output.failed}
cat {input.sup} | awk -v asm="{wildcards.asm}" -v meth="{wildcards.meth}" -v hap="{wildcards.hap}" '{{if ($7 >= 3 || $8 >= 3) {{print $1,$2,$3,asm"."hap"/"meth"/"$4,$5,$6,$7,$8 }} }}' OFS="\t" > tmp.{wildcards.asm}.{wildcards.meth}.{wildcards.hap}.passed

cat tmp.{wildcards.asm}.{wildcards.meth}.{wildcards.hap}.passed | awk -v avg="$avgCov"\
                      '{{ cmd="samtools coverage -r "$1":"$2"-"$3" {input.reads} | tail -n1 | cut -f7"; cmd | getline localCov ; if(localCov < 2*avg) {{print}} }}' OFS="\t" > {output.passed} && rm -f tmp.{wildcards.asm}.{wildcards.meth}.{wildcards.hap}.passed
"""

#
#Step 8: Collect statistics on this run of the pipeline
#
rule CollectStats:
    input:
        asm="indices/{asm}.{hap}.assembly.fasta.fai",
        covs=expand("{{asm}}.{{hap}}/{method}/aln.bam", method=methods),
        vars=expand("{{asm}}.{{hap}}/{method}/variants.sv.bed", method=methods),
        sups=expand("{{asm}}.{{hap}}/{method}/sv.passed.bed", method=methods) 
    output:
        stats="{asm}.{hap}/stats.txt",
    params:
        grid_opts=config["grid_small"],
        ref=config["ref"]
    shell:"""
covFiles[0]={input.covs[0]};covFiles[1]={input.covs[1]} ; varFiles[0]={input.vars[0]};varFiles[1]={input.vars[1]} ; supFiles[0]={input.sups[0]};supFiles[1]={input.sups[1]}

halfGnm=$(cat {params.ref}.fai | awk '{{size += $2}} END {{print int(size/2)}}')

n50=$(sort -nr -k 2 {input.asm} | awk -v gnm="$halfGnm" '{{contgSum += $2 ; if(contgSum >= gnm) {{print $2}} }} END{{if(contgSum < gnm) {{print "notFound"}} }}' \ 
    > {output.stats}.tmp ; head -n 1 {output.stats}.tmp) && rm -f {output.stats}.tmp
for i in {{0..1}} ; do covs[$i]=$(bedtools genomecov -ibam ${{covFiles[$i]}} | awk '{{if ($1=="genome" && $2>0) {{count += $5}} }} END {{print count}}') ; done
for i in {{0..1}} ; do vars[$i]=$(cat ${{varFiles[$i]}} | wc -l) ; done
for i in {{0..1}} ; do sups[$i]=$(cat ${{supFiles[$i]}} | wc -l) ; done
echo -e "Assembly\tN50\tmm2BasesCov\tlraBasesCov\tmm2SVs\tlraSVs\tmm2Sup\tlraSup" > {output.stats}
echo -e "{wildcards.asm}.{wildcards.hap}\t$n50\t${{covs[1]}}\t${{covs[0]}}\t${{vars[1]}}\t${{vars[0]}}\t${{sups[1]}}\t${{sups[0]}}" >> {output.stats}
"""

rule CombineStats:
    input:
        stats=expand("{asm}.{hap}/stats.txt",asm=asms,hap=haps)
    output:
        all="run_stats_all.txt"
    shell:"""

echo -e "Assembly\tN50\tmm2BrdthCov\tlraBrdthCov\tmm2SVs\tlraSVs\tmm2Sup\tlraSup" > {output.all}
for inpStat in {input.stats}; do cat $inpStat | tail -n 1 >> {output.all} ; done
"""
