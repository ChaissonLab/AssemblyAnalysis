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
    asmFiles[asm]["1"] = config["assemblies"][asm][0]
    asmFiles[asm]["2"] = config["assemblies"][asm][1]

methods=["lra", "mm2"]    
ops=["ins", "del"]
cats=["missed", "found"]
cats +=["missed_slop", "found_slop"]
samples=config["sample"]

localrules:CollectStats,CompareToSVSet,GetVariantSupport,CombineVariantSupport,ChromVCFToBed,CombineOps,CombOps,SplitFasta,LraBed,LraClust,LraClustDip,SplitLrabed

rule all:
    input:
        gapVCF=expand("from_reads.{method}.vcf",method=methods),
        gapBeds=expand("from_reads.{method}.{op}.bed",method=methods,op=ops),
        lrasplit=expand("{asm}.{hap}/lra/split.fasta", asm=asms,hap=haps),
        lrabed=expand("{asm}.{hap}/lra/calls.bed", asm=asms, hap=haps),
        lrabedclust=expand("{asm}.{hap}/lra/calls.clust.bed", asm=asms, hap=haps),
        lraDipClust=expand("{asm}_dip/lra/variants.clustered.bed",asm=asms),
        lrabedop=expand("{asm}.{hap}/lra/variants.sv.{op}.bed", asm=asms, hap=haps, op=ops),
        lrabedopSup=expand("{asm}.{hap}/lra/variants.sv.{op}.bed.{meth}.support", meth=methods, asm=asms, hap=haps, op=ops),
        lrabedopSupComb=expand("{asm}.{hap}/lra/variants.sv.{op}.bed.support.combined", asm=asms, hap=haps, op=ops),
        sam=expand("{asm}.{hap}/lra/split.fasta.bam", asm=asms,hap=haps),
        mm2HapVcf=expand("{asm}.{hap}/mm2/variants.vcf", asm=asms, hap=haps),
        paf2bed=expand("{asm}.{hap}/mm2/aln.paf.bed", asm=asms,hap=haps),
        mm2HapVcfGz=expand("{asm}.{hap}/mm2/variants.vcf.bgz", asm=asms, hap=haps),
        mm2HapSVVcf=expand("{asm}.{hap}/mm2/variants.sv.vcf", asm=asms, hap=haps),
        mm2HapSVBedAll=expand("{asm}.{hap}/mm2/variants.sv.bed", asm=asms, hap=haps),
        mm2HapSVBed=expand("{asm}.{hap}/mm2/variants.sv.{op}.bed", asm=asms, hap=haps,op=ops),
        mm2HapSVBedSup=expand("{asm}.{hap}/mm2/variants.sv.{op}.bed.{method}.support", method=methods,asm=asms, hap=haps,op=ops),
        mm2HapSVBedSupComb=expand("{asm}.{hap}/mm2/variants.sv.{op}.bed.support.combined", asm=asms, hap=haps,op=ops),
        mm2HapVcfToDip=expand("{asm}_dip/mm2/variants.vcf.gz", asm=asms,op=ops),
        mm2Comb=expand("{asm}.{hap}/mm2/variants.clusters.bed", asm=asms, hap=haps),
        mm2CombDip=expand("{asm}_dip/mm2/variants.clustered.bed", asm=asms, hap=haps),
        comparison=expand("{asm}_dip/{meth}/variants.sv.{op}.{cat}.bed", asm=asms, op=ops+["hap"], cat=cats, meth=methods),
        stats=expand("run_stats.txt", asm=asms, hap=haps, meth=methods)

#        mm2DipSVBed=expand("{asm}_dip/mm2/variants.sv.{op}.bed", asm=asms,op=ops),



rule GetVariantSupport:
    input:
        bed="{prefix}.{op}.bed",
        sup="from_reads.{method}.{op}.bed"
    output:
        sup="{prefix}.{op}.bed.{method}.support"
    resources:
        threads=1
    params:
        ref=config["ref"],
        sd=SD
    shell:"""
bedtools slop -g {params.ref}.fai -b 1000 -i {input.sup} | \
  bedtools intersect -loj -a <( cut -f 1-5 {input.bed} ) -b stdin | \
  bedtools groupby -c 10 -o collapse -full | \
  {params.sd}/GetSVSupport.py > {output.sup}
"""

rule CombineVariantSupport:
    input:
        meth=expand("{{prefix}}.{{op}}.bed.{method}.support", method=methods),
    output:
        comb="{prefix}.{op}.bed.support.combined"
    shell:"""
paste {input.meth[0]} <( cut -f 6 {input.meth[1]} ) > {output.comb}
"""

rule CreateChromVCF:
    input:
        reads=lambda wildcards: config["bam"][wildcards.method]
    output:
        vcf="gaps/{method}.{chrom}.gaps.vcf"
    resources:
        threads=16
    params:
        grid_opts=config["grid_manycore"],
        ref=config["ref"],
        sample=config["sample"],
        sd=SD
    shell:"""
mkdir -p gaps
{params.sd}/SamToVCF.py --ref {params.ref} --sample {params.sample} --sam {input.reads} --minLength 50 --chrom {wildcards.chrom} > {output.vcf}
"""

rule CreateGapBed:
    input:
        vcfs=expand("gaps/{{method}}.{chrom}.gaps.vcf", chrom=allChroms)
    output:
        gapVcf="from_reads.{method}.vcf"
    resources:
        threads=16
    params:
        grid_opts=config["grid_manycore"],
        ref=config["ref"],
        sample=config["sample"],
        sd=SD
    shell:"""
grep "^#" {input.vcfs[0]} > {output.gapVcf}
cat {input.vcfs} | grep -v "^#" >> {output.gapVcf}
"""

rule ChromVCFToBed:
    input:
        vcf="from_reads.{method}.vcf"
    output:
        bed=expand("from_reads.{{method}}.{op}.bed",op=ops)
    params:
        sd=SD
    shell:"""
{params.sd}/VcfToSplitBed.sh {input.vcf} {output.bed[1]} {output.bed[0]}
"""

rule CombineOps:
    input:
        opBed="{asm}.{hap}/mm2/variants.sv.bed"
    output:
        bed="{asm}.{hap}/mm2/variants.clusters.bed"
    shell:"""
cat {input.opBed} | bedtools cluster -d 500 | bedtools groupby -g 7 -c 1,2,3,4,5,6,7 -o first,min,max,collapse,collapse,collapse,count | cut -f 2- | \
awk '{{ if (index($6,"deletion") > 0 && index($6,"insertion") > 0) {{ $6="locus";}} else {{ if (index($6,"deletion") == 0) {{ $6="insertion"; }} else {{ $6="deletion";}} }} print; }}' | tr " " "\\t"  > {output.bed}
"""

rule CombOps:
    input:
        bed=expand("{{asm}}.{hap}/mm2/variants.clusters.bed",hap=haps)
    output:
        mm2CombDip="{asm}_dip/mm2/variants.clustered.bed",
    shell:"""
cat {input.bed} | bedtools sort | bedtools cluster -d 500 | bedtools groupby -c 8 -g 8 -o count -full > {output.mm2CombDip}
"""
rule SplitFasta:
    input:
        asm=lambda wildcards: asmFiles[wildcards.asm][wildcards.hap]
    output:
        split="{asm}.{hap}/lra/split.fasta",
        index="{asm}.{hap}/assembly.fasta.fai"
    params:
        splitSize=500000,
        sd=SD
    shell:"""
ln -s {input.asm} {wildcards.asm}.{wildcards.hap}/assembly.fasta && samtools faidx {wildcards.asm}.{wildcards.hap}/assembly.fasta && rm {wildcards.asm}.{wildcards.hap}/assembly.fasta

mkdir -p {wildcards.asm}.{wildcards.hap}/lra
{params.sd}/SplitFasta.py {input.asm} {params.splitSize} > {output.split}
"""

rule MapSplitFasta:
    input:
        split="{asm}.{hap}/lra/split.fasta"
    output:
        sam="{asm}.{hap}/lra/split.fasta.bam"
    params:
        ref=config["ref"],
        grid_opts=config["grid_large"]
    shell:"""
lra align {params.ref} {input.split} -t 8 -p s | \
 samtools sort -@2 -T $TMPDIR/$$  -o {output.sam}
samtools index {output.sam}

"""
rule LraBed:
    input:
        bam="{asm}.{hap}/lra/split.fasta.bam"
    output:
        bed="{asm}.{hap}/lra/calls.bed"
    params:
        ref=config["ref"],
        sd=SD
    shell:"""
samtools view -h {input.bam} | {params.sd}/PrintGaps.py {params.ref} /dev/stdin | awk -F " " '!keep[$1,$2,$3,$4]++' > {output.bed}
"""

rule LraClust:
    input:
        bed="{asm}.{hap}/lra/calls.bed"
    output:
        lrabed="{asm}.{hap}/lra/calls.clust.bed"
    shell:"""
bedtools sort -header -i {input.bed} | \
  bedtools cluster -d 500 | \
  bedtools groupby -g 11 -c 1,2,3,4,5,6,11 -o first,min,max,collapse,collapse,collapse,count | \
  cut -f 2- | \
  awk '{{ if (index($4,"deletion") > 0 && index($4,"insertion") > 0) {{ $4="locus";}} else {{ if (index($4,"deletion") == 0) {{ $4="insertion"; }} else {{ $4="deletion";}} }} print; }}' | \
  tr " " "\\t"  > {output.lrabed}
"""

rule LraClustDip:
    input:
        lrabed=expand("{{asm}}.{hap}/lra/calls.clust.bed",hap=haps)
    output:
        dip="{asm}_dip/lra/variants.clustered.bed"
    shell:"""
cat {input.lrabed} | bedtools sort | bedtools cluster -d 500 | bedtools groupby -c 8 -g 8 -o count -full > {output.dip}
"""

rule SplitLrabed:
    input:
        bed="{asm}.{hap}/lra/calls.bed",
    output:
        split="{asm}.{hap}/lra/variants.sv.{op}.bed",
    shell:"""

egrep "^#|insertion" {input.bed} > {wildcards.asm}.{wildcards.hap}/lra/variants.sv.ins.bed
egrep "^#|deletion" {input.bed} > {wildcards.asm}.{wildcards.hap}/lra/variants.sv.del.bed
"""
rule CompareToSVSet:
    input:
        bed="{asm}_dip/{method}/variants.clustered.bed",
        other=expand("{samp}.hgsvg-p1.{{op}}.bed",samp=samples)
    output:
        missed="{asm}_dip/{method}/variants.sv.{op}.missed.bed",
        found="{asm}_dip/{method}/variants.sv.{op}.found.bed",
        missed_slop="{asm}_dip/{method}/variants.sv.{op}.missed_slop.bed",
        found_slop="{asm}_dip/{method}/variants.sv.{op}.found_slop.bed"

    params:
        grid_opts=config["grid_blat"],
        ref=config["ref"]
    shell:"""

bedtools intersect -v -a {input.other} -b {input.bed} > {output.missed}
bedtools intersect -u -a {input.other} -b {input.bed} > {output.found}

bedtools slop -i {input.bed} -b 1000 -g {params.ref}.fai | bedtools intersect -v -a {input.other} -b stdin > {output.missed_slop}
bedtools slop -i {input.bed} -b 1000 -g {params.ref}.fai | bedtools intersect -u -a {input.other} -b stdin > {output.found_slop}
"""

rule PafToBed:
    input:
        paf="{asm}.{hap}/mm2/aln.paf"
    output: 
        bed="{asm}.{hap}/mm2/aln.paf.bed"
    params:
        grid_opts=config["grid_small"]
    shell:"""
cut -f 6,8,9 {input.paf} | bedtools sort | bedtools merge >  {output.bed}
"""

rule hapvcfToDip:
    input:
        vcf=expand("{{asm}}.{hap}/{{method}}/variants.vcf.bgz", hap=haps)
    output:
        dipvcf="{asm}_dip/{method}/variants.vcf.gz"
    params:
        grid_opts=config["grid_small"],    
    shell:"""
mkdir -p {wildcards.asm}_dip/{wildcards.method}

bcftools merge -0 --force-samples {input.vcf}  | awk 'BEGIN {{OFS="\\t"; }} {{ if (substr($0,0,1) == "#") {{ print; }} else {{ if ($(NF-1) != $NF) {{ $(NF-1) = "0/1";}} print;}} }}' | cut -f 1-10 | bgzip -c > {output.dipvcf}


"""
    
    
rule bgzipvcf:
    input:
        vcf="{asm}.{hap}/{method}/variants.vcf"
    output:
        bgvcf="{asm}.{hap}/{method}/variants.vcf.bgz"
    params:
        grid_opts=config["grid_small"],    
    shell:"""
bgzip -c {input.vcf} > {output.bgvcf}
tabix {output.bgvcf}
"""




# 
# First step in minimap2 pipeline: map contigs, and call variants using Heng's paftools pipeline
# 
rule MakeMM2VCF:
    input:
        asm=lambda wildcards: asmFiles[wildcards.asm][wildcards.hap]
    output:
        vcf="{asm}.{hap}/mm2/variants.vcf",
        paf="{asm}.{hap}/mm2/aln.paf"
    params:
        grid_opts=config["grid_large"],
        ref=config["ref"]
    shell:"""


mkdir -p {wildcards.asm}.{wildcards.hap}/mm2
minimap2 -t 16 -x asm5 --cs {params.ref} {input.asm}  > {output.paf}
cat {output.paf} | sort -k6,6 -k8,8n | \
  paftools.js call -f {params.ref} - > {output.vcf}

"""

#
# Step 2: generate SV calls from full minimap2 pipeline -- calls at least 50 bases.
#

rule SVVCF:
    input:
        vcf="{asm}.{hap}/mm2/variants.vcf"
    output:
        sv="{asm}.{hap}/mm2/variants.sv.vcf"
    params:
        grid_opts=config["grid_small"],
    shell:"""
cat {input.vcf} | awk '{{ if (substr($1,0,1) == "#") {{ print ; }} else {{ a=length($4); b=length($5); diff=a-b; if (diff <= -50 || diff >= 50) print;}} }}' > {output.sv}
"""

#
# Step 3: Create bed files with SV calls split by operation: deletion and insertion. 
# ************* CAUTION THIS SHOULD BE IMPROVED ****************
#

rule SVVCFToBed:
    input:
        vcf="{asm}.{hap}/mm2/variants.sv.vcf"
    output:
        vcf=expand("{{asm}}.{{hap}}/mm2/variants.sv.{op}.bed",op=ops),
        bed="{asm}.{hap}/mm2/variants.sv.bed"
    params:
        grid_opts=config["grid_small"],
    shell:"""


grep -v "^#" {input.vcf} | awk '{{ a=length($4); b=length($5); diff=a-b; {{if (diff <= 0) {{oper="insertion"; diff*= -1}} else {{oper="deletion"}};\
                                 print $1,$2,$2+diff,$4,$5,oper; !keep[$1,$2,$3,$6]++}} }}' OFS="\t" | bedtools sort > {output.bed} 

cat {output.bed} | awk '{{if ($6 == "insertion") {{print}}  }}'  >  {wildcards.asm}.{wildcards.hap}/mm2/variants.sv.ins.bed
cat {output.bed} | awk '{{if ($6 == "deletion")  {{print}}  }}'  >  {wildcards.asm}.{wildcards.hap}/mm2/variants.sv.del.bed
"""

rule CollectStats:
    input:
         asm=expand("{asm}.{hap}/assembly.fasta.fai",asm=asms,hap=haps),
         covL=expand("{asm}.{hap}/lra/split.fasta.bam",asm=asms,hap=haps),
         covM=expand("{asm}.{hap}/mm2/aln.paf",asm=asms,hap=haps),
         varL=expand("{asm}.{hap}/lra/calls.bed",asm=asms,hap=haps),
         varM=expand("{asm}.{hap}/mm2/variants.sv.bed",asm=asms,hap=haps),
         sup=expand("{asm}.{hap}/{meth}/variants.sv.{op}.bed.support.combined",asm=asms,hap=haps,meth=methods,op=ops)
    output:
        stats="run_stats.txt",
        false="falsePosAll.bed"
    params:
        ref=config["ref"]
    shell:"""

echo -e "Assembly\tmm2SVs\tlraSVs" > {output.stats}
for inpM in {input.varM}; do wc -l $inpM | awk '{{split($2, arr, "/") ; assemb=arr[1] ; print assemb,$1,"."}}' OFS="\t" >> {output.stats}  ; done
for inpL in {input.varL}; do svCount=$(cat $inpL | wc -l) ; asmName=$(echo $inpL | cut -d'/' -f1) ;\
                          cat {output.stats} | awk -v count="$svCount" -v name="$asmName" '{{if ($1 == name) {{ print $1,$2,count}} else {{print $0}}  }}' OFS="\t" \
                          > {output.stats}.tmp && mv -f {output.stats}.tmp {output.stats} ; done
cp {output.stats} {output.stats}.tmp ; echo "" > {output.false}
for inpS in {input.sup};  do supCount=$(cat $inpS | awk '{{if ($6 >= 3 || $7 >= 3) {{print}} }}' | wc -l) ; asmName=$(echo $inpS | cut -d'/' -f1) ;\
                          methName=$(echo $inpS | cut -d'/' -f2) ; opName=$(echo $inpS | cut -d'/' -f3 | cut -d'.' -f3) ;\
                          cat $inpS | awk -v name="$asmName" -v opr="$opName" -v met="$methName" '{{if ($6 < 3 && $7 < 3) {{print $1,$2,$3,name"/"met"/"opr}} }}' OFS="\t" >> {output.false} ;\
                          cat {output.stats}.tmp | awk -v count="$supCount" -v name="$asmName" -v met="$methName" -v opr="$opName"\
                          '{{if($1 == name)    {{if (met == "mm2" && opr == "ins") {{print $1,count,$3,$4,$5}}\
                                                 else if(met == "mm2" && opr == "del") {{print $1,$2,count,$4,$5}}\
                                                 else if(met == "lra" && opr == "ins") {{print $1,$2,$3,count,$5}}\
                                                 else if(met == "lra" && opr == "del") {{print $1,$2,$3,$4,count}} }} \
                                                 else {{print}} }}' OFS="\t" \
                                                 > {output.stats}.tmp.tmp && mv -f {output.stats}.tmp.tmp {output.stats}.tmp ; done
cat {output.stats}.tmp | awk 'BEGIN{{print "mm2Sup\tlraSup"}} NR>1 {{print $2+$3,$4+$5}}' OFS="\t" | paste {output.stats} - \
                          > {output.stats}.tmp.tmp && mv -f {output.stats}.tmp.tmp {output.stats} && rm -f {output.stats}.tmp 
cp {output.stats} {output.stats}.tmp
for inpM in {input.covM}; do breadth=$(bedtools genomecov -i $inpM.bed -g {params.ref}.fai | awk '{{if ($1=="genome" && $2>0) {{count += $5}}}} END {{print count}}') ; asmName=$(echo $inpM | cut -d'/' -f1) ;\
                          cat {output.stats}.tmp | awk -v brdth="$breadth" -v name="$asmName" \
                          '{{if($1 == name) {{print $1, brdth}} else {{print}} }}' OFS="\t"\
                          > {output.stats}.tmp.tmp && mv -f {output.stats}.tmp.tmp {output.stats}.tmp ; done
cat {output.stats}.tmp | awk 'BEGIN{{print "mm2BrdthCov"}} NR>1 {{print $2}}' OFS="\t" | paste {output.stats} - \
                          > {output.stats}.tmp.tmp && mv -f {output.stats}.tmp.tmp {output.stats} && rm -f {output.stats}.tmp
cp {output.stats} {output.stats}.tmp
for inpL in {input.covL}; do breadth=$(bedtools genomecov -ibam $inpL | awk '{{if ($1=="genome" && $2>0) {{count += $5}}}} END {{print count}}') ; asmName=$(echo $inpL | cut -d'/' -f1) ;\
                          cat {output.stats}.tmp | awk -v brdth="$breadth" -v name="$asmName" \
                          '{{if($1 == name) {{print $1, brdth}} else {{print}} }}' OFS="\t"\
                          > {output.stats}.tmp.tmp && mv -f {output.stats}.tmp.tmp {output.stats}.tmp ; done
cat {output.stats}.tmp | awk 'BEGIN{{print "lraBrdthCov"}} NR>1 {{print $2}}' OFS="\t" | paste {output.stats} - \
                          > {output.stats}.tmp.tmp && mv -f {output.stats}.tmp.tmp {output.stats} && rm -f {output.stats}.tmp

echo "N50" > {output.stats}.tmp.tmp
cp {output.stats} {output.stats}.tmp
for assemb in {input.asm}; do gnmSize=$(cat {params.ref}.fai | awk '{{size += $2}} END {{print int(size/2)}}') ; sort -nr -k 2 $assemb | awk -v gnm="$gnmSize" \
                          'BEGIN{{contgSum=0}} {{contgSum += $2 ; if(contgSum >= gnm) {{print $2}} }}' > tmp ; head -n 1 tmp >> {output.stats}.tmp.tmp ; done ; paste {output.stats}.tmp {output.stats}.tmp.tmp \
                          > {output.stats} ; rm -f {output.stats}.tmp {output.stats}.tmp.tmp tmp


"""
