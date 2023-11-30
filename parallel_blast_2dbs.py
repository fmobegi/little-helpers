#!/usr/bin/env python3
# File created on 06 Feb 2023

"""
This program takes two blast databases and a query sequence, 
performs blastn seraches and only return the best matched database if it has more than 10% of the total matches.
BLAST output results are filtered for a minimum sequence identity of 95% 
and a minimum subject coverage of 80% of the subject sequence. 
"""
import argparse
import multiprocessing
import os
import subprocess
from Bio import SeqIO

def run_blastn(query_file, output_file, database):
    blast_output = subprocess.run(['blastn', '-query', query_file, '-db', database, '-outfmt', '6 qseqid sseqid pident bitscore evalue length mismatch staxids sscinames scomnames sskingdoms stitle', '-max_target_seqs', '1', '-max_hsps', '1', '-evalue', '1e-10'], stdout=subprocess.PIPE).stdout.decode()
    lines = blast_output.strip().split('\n')
    filtered_output = []
    for line in lines:
        fields = line.split('\t')
        pident = float(fields[2])
        length = int(fields[3])
        qlen = int(fields[12])
        if pident >= 90.0 and length/qlen >= 0.8:
            filtered_output.append(line)
    with open(output_file, 'w') as f:
        f.write('\n'.join(filtered_output))

def compare_results(results):
    # Compare the results of both BLAST runs and return the best match
    result1, result2 = results
    best_match = None
    count1 = sum(1 for r in result1 if r[1] > 0)
    count2 = sum(1 for r in result2 if r[1] > 0)
    total_count = count1 + count2
    if total_count == 0:
        return None
    if count1/total_count >= 0.1:
        best_match = "database1"
    elif count2/total_count >= 0.1:
        best_match = "database2"
    return best_match

def main(args):
    # Read the query sequence from the FASTA file
    records = SeqIO.parse(args.input, "fasta")
    query_records = list(records)

    results = []
    # Start a process pool with number of workers equal to the number of cores
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        for _ in query_records:
            basename, _ = (
                os.path.splitext(args.output)
                if args.output
                else os.path.splitext(args.input)
            )
            output_file1 = f'{basename}_db1.out'
            output_file2 = f'{basename}_db2.out'
            results.extend(
                (
                    pool.apply_async(
                        run_blastn,
                        args=(args.input, output_file1, args.database1),
                    ),
                    pool.apply_async(
                        run_blastn,
                        args=(args.input, output_file2, args.database2),
                    ),
                )
            )
        pool.close()
        pool.join()

    best_match = compare_results(results)
    if best_match is None:
        print("No database returned matches for at least 10 percent of the query sequences.")
    else:
        print("Best matching database:", best_match)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parallel BLASTn")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file name")
    parser.add_argument("-d1", "--database1", required=True, help="BLAST database name 1")
    parser.add_argument("-d2", "--database2", required=True, help="BLAST database name 2")
    parser.add_argument("-o", "--output", help="Output file name (default: file basename + .out)")
    args = parser.parse_args()
    main(args)
