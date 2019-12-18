import os
import tempfile
import subprocess
import os.path
env = os.environ
import sys

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

# Snakemake and working directories
SD = os.path.dirname(workflow.snakefile)

# Config
configfile: "multi_asm.json"

ref=config["ref"]

fai = open(ref + ".fai")
allChroms = [l.strip().split()[0] for l in fai ]

haps=["1", "2"]
assemblers=config["assemblers"]

datasets = [entry for entry in config["datasets"].keys()]

rule all:
    input:
        config=expand("{dataset}/{assembler}/partitioned_assembly.json", dataset=datasets, assembler=assemblers, hap=haps),
        asm=expand("{dataset}/{assembler}/assembly.{hap}.consensus.fasta", dataset=datasets, assembler=assemblers, hap=haps)


rule GenerateJSON:
    input:
        bams=lambda wildcards: config["datasets"][wildcards.datasetID]["bams"],
        vcf=lambda wildcards: config["datasets"][wildcards.datasetID]["vcf"]
    output:
        json="{datasetID}/{assembler}/partitioned_assembly.json",
    params:
        sd=SD,
        sample=lambda wildcards: config["datasets"][wildcards.datasetID]["sample"],
        wd=lambda wildcards: config["workingDir"] + "/" + wildcards.datasetID,
        ref=config["ref"],
        readtype=lambda wildcards: config["datasets"][wildcards.datasetID]["datatype"],
        workingDir=config["workingDir"]
    shell:"""
mkdir -p {wildcards.datasetID}/{wildcards.assembler}
{params.sd}/CreateAssemblyJSON.py --bams {input.bams} --workingDir {params.workingDir}/{wildcards.datasetID}/{wildcards.assembler}/run --sample {params.sample} --ref {params.ref} --vcf {input.vcf} --readtype {params.readtype} --assembler {wildcards.assembler} > {wildcards.datasetID}/{wildcards.assembler}/partitioned_assembly.json
"""


rule RunAssembly:
    input:
        json="{dataset}/{assembler}/partitioned_assembly.json",
    output:
        asms=expand("{{dataset}}/{{assembler}}/assembly.{hap}.consensus.fasta",hap=haps)
    params:
        wd=lambda wildcards: config["workingDir"] + "/" + wildcards.dataset  + "/" + wildcards.assembler,
        pd=config["pd"],
        jobs_per_run=config["jobs_per_run"]
    shell:"""
mkdir -p {params.wd}
cp {input.json} {params.wd}/
pushd {params.wd} && snakemake -p -s {params.pd}/PartitionedAssembly.snakefile  -j {params.jobs_per_run} --cluster " {{params.grid_opts}} -c {{resources.threads}} --mem={{resources.mem_gb}}G {{params.node_constraint}} "  --restart-times 4  && popd
cp {params.wd}/assembly.*.consensus.fasta {wildcards.dataset}/{wildcards.assembler}/
rm -rf {params.wd}/

"""
