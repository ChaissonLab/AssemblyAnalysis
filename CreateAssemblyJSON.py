#!/usr/bin/env python
import sys
import argparse
ap=argparse.ArgumentParser(description="Create a .json file for running an assembly.")
ap.add_argument("--assembler", help="Which assembler to use: wtdbg2, flye, canu, or falcon", required=True)
ap.add_argument("--readtype", help="Read type: pacbio, ccs, or ont", required=True)
ap.add_argument("--bams", help="Input bams", nargs="+", required=True)
ap.add_argument("--workingDir", help="Run the assembly in this directory", required=True)
ap.add_argument("--sample", help="Sample", required=True)
ap.add_argument("--vcf", help="Phasing vcf", required=True)
ap.add_argument("--ref", help="Ref", required=True)
ap.add_argument("--partition", help="Optional partition", default="cmb")
args=ap.parse_args()

if args.readtype == "ccs" or args.readtype == "ont":
    consensus = "racon"
else:
    consensus = "arrow"

part=args.partition

sys.stdout.write("""{{ 
    "bam" : [{}],
    "working_dir" : "{}",
    "sample" : "{}",
    "ref" : "{}",
    "vcf" : "{}",
    "partition" : {},
    "grid_large" : "sbatch --time=24:00:00 --partition={}",
    "grid_manycore" : "sbatch --time=48:00:00 --partition={}",
    "grid_medium" : "sbatch --time=24:00:00 --partition={}",
    "grid_small" : "sbatch --time=4:00:00 --partition={}",
    "grid_blat" : "sbatch --time=24:00:00 --partition={}",
    "read-type" : "{}",
    "consensus" : "{}",
    "assembler" : "{}"
}}""".format(", ".join(args.bams), args.workingDir, args.sample, args.ref, args.vcf, 
             args.partition,
             args.partition,
             args.partition,
             args.partition,
             args.partition,
             args.partition, args.readtype, consensus, args.assembler))

