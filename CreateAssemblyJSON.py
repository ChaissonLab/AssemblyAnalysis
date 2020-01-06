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
args=ap.parse_args()

if args.readtype == "ccs" or args.readtype == "ont":
    consensus = "racon"
else:
    consensus = "arrow"

sys.stdout.write("""{{ 
    "bam" : [{}],
    "working_dir" : "{}",
    "sample" : "{}",
    "ref" : "{}",
    "vcf" : "{}",
    "grid_large" : "sbatch --time=24:00:00 --partition=cmb",
    "grid_manycore" : "sbatch --time=48:00:00 --partition=cmb",
    "grid_medium" : "sbatch --time=24:00:00 --partition=cmb",
    "grid_small" : "sbatch --time=4:00:00 --partition=cmb",
    "grid_blat" : "sbatch --time=24:00:00 --partition=cmb",
    "read-type" : "{}",
    "consensus" : "{}",
    "assembler" : "{}"
}}""".format(", ".join(args.bams), args.workingDir, args.sample, args.ref, args.vcf, args.readtype, consensus, args.assembler))

