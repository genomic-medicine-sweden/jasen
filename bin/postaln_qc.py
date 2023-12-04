import argparse
import sys
import os
import subprocess
import json

description = '''
Run QC after alignment.
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog="python3 postaln_qc.py --bam ${bam} --ref_fasta ${reference} --sample_id ${sampleName} --threads ${task.cpus} --output ${output}"
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.0.1'
    )
parser.add_argument(
    "--bam",
    help="Input BAM file.",
    dest='bam',
    required=True
    )
parser.add_argument(
    "--bed",
    help="Targeted BED file.",
    dest='bed'
    )
parser.add_argument(
    "--sample_id",
    help="Sample ID.",
    dest='sample_id',
    required=True
    )
parser.add_argument(
    "--threads",
    help="Number of threads.",
    type=int,
    default=0,
    dest='threads'
    )
parser.add_argument(
    "--baits",
    help="Baits BED file.",
    default=0,
    dest='baits'
    )
parser.add_argument(
    "--ref_fasta",
    help="Reference FASTA file.",
    default=0,
    dest='ref_fasta',
    required=True
    )
parser.add_argument(
    '-o', '--output',
    help='output filepath',
    metavar='OUTPUT_FPATH',
    dest='output',
    required=True
    )

args = parser.parse_args()

class QC:
    def __init__(self, args):
        self.results = {}
        self.bam = args.bam
        self.bed = args.bed
        self.sample_id = args.sample_id
        self.threads = args.threads
        self.baits = args.baits
        self.ref_fasta = args.ref_fasta
        self.paired = self.is_paired()

    def run_qc(self):
        if self.baits and self.ref_fasta:
            print("Calculating HS-metrics...")
            dict_file = self.ref_fasta
            if not dict_file.endswith(".dict"):
                dict_file += ".dict"
            if not os.path.isfile(f"{self.bed}.interval_list"):
                self.system_p(f"picard BedToIntervalList -I {self.bed} -O {self.bed}.interval_list -SD {dict_file}")
            if not os.path.isfile(f"{self.baits}.interval_list"):
                self.system_p(f"picard BedToIntervalList -I {self.baits} -O {self.baits}.interval_list -SD {dict_file}")
            self.system_p(f"picard CollectHsMetrics -I {self.bam} -O {self.bam}.hsmetrics -R {self.ref_fasta} -BAIT_INTERVALS {self.baits}.interval_list -TARGET_INTERVALS {self.bed}.interval_list")

            with open(f"{self.bam}.hsmetrics") as hs:
                for line in hs:
                    if line.startswith("## METRICS CLASS"):
                        next(hs)
                        vals = next(hs).split("\t")
                        self.results['pct_on_target'] = vals[18]
                        self.results['fold_enrichment'] = vals[26]
                        self.results['median_coverage'] = vals[23]
                        self.results['fold_80'] = vals[33]

        print("Collecting basic stats...")
        flagstat = subprocess.check_output(f"sambamba flagstat {'-t '+str(self.threads) if self.threads else ''} {self.bam}", shell=True, text=True).splitlines()
        num_reads = int(flagstat[0].split()[0])
        dup_reads = int(flagstat[3].split()[0])
        mapped_reads = int(flagstat[4].split()[0])

        if self.paired:
            print("Collect insert sizes...")
            self.system_p(f"picard CollectInsertSizeMetrics -I {self.bam} -O {self.bam}.inssize -H {self.bam}.ins.pdf -STOP_AFTER 1000000")
            with open(f"{self.bam}.inssize") as ins:
                for line in ins:
                    if line.startswith("## METRICS CLASS"):
                        next(ins)
                        vals = next(ins).split("\t")
                        self.results['ins_size'] = vals[0]
                        self.results['ins_size_dev'] = vals[1]

            os.remove(f"{self.bam}.inssize")
            os.remove(f"{self.bam}.ins.pdf")

        out_prefix = f"{self.bam}_postalnQC"
        thresholds = [1, 10, 30, 100, 250, 500, 1000]

        print("Collecting depth stats...")
        self.system_p(f"sambamba depth base -c 0 {'-t '+str(self.threads) if self.threads else ''} -L {self.bed} {self.bam} > {out_prefix}.basecov.bed")
        pct_above, mean_cov, iqr_median = self.parse_basecov_bed(f"{out_prefix}.basecov.bed", thresholds)
        os.remove(f"{out_prefix}.basecov.bed")

        self.results['pct_above_x'] = pct_above
        self.results['tot_reads'] = num_reads
        self.results['mapped_reads'] = mapped_reads
        self.results['dup_reads'] = dup_reads
        self.results['dup_pct'] = dup_reads / mapped_reads
        self.results['sample_id'] = self.sample_id
        self.results['mean_cov'] = mean_cov
        self.results['iqr_median'] = iqr_median

        json_result = json.dumps(self.results, indent=4)
        return json_result

    def write_json_result(self, json_result, output_filepath):
        with open(output_filepath, 'w') as json_file:
            json_file.write(json_result)

    def parse_basecov_bed(self, fn, thresholds):
        with open(fn) as cov_fh:
            head_str = cov_fh.readline().strip().lstrip("#")
            head = head_str.split("\t")
            cov_field = head.index("COV")

            tot_bases = 0
            above_cnt = {min_val: 0 for min_val in thresholds}

            tot, cnt = 0, 0
            levels = {}
            for line in cov_fh:
                a = line.strip().split("\t")
                tot += int(a[2])
                cnt += 1
                tot_bases += 1
                for min_val in thresholds:
                    if int(a[cov_field]) >= min_val:
                        above_cnt[min_val] += 1

            above_pct = {min_val: 100 * (above_cnt[min_val] / tot_bases) for min_val in thresholds}

            mean_cov = tot / cnt

            # Calculate the inter-quartile range / median (IQR/median)
            q1_num = cnt / 4
            q3_num = 3 * cnt / 4
            median_num = cnt / 2
            sum_val = 0
            q1, q3, median = None, None, None
            iqr_median = "9999"
            for l in sorted(levels):
                sum_val += levels[l]
                if sum_val >= q1_num and not q1:
                    q1 = l
                if sum_val >= median_num and not median:
                    median = l
                if sum_val >= q3_num and not q3:
                    q3 = l

            if q1 and q3 and median:
                iqr_median = (q3 - q1) / median

            return above_pct, mean_cov, iqr_median

    def is_paired(self):
        line = subprocess.check_output(f"samtools view {self.bam} | head -n 1| awk '{{print $2}}'", shell=True, text=True)
        remainder = int(line) % 2
        is_paired = 1 if remainder else 0
        return is_paired

    def system_p(self, *cmd):
        print(f"RUNNING: {' '.join(cmd)}")
        print()
        subprocess.run(cmd, check=True)

def main():
    qc = QC(args)
    json_result = qc.run_qc()
    qc.write_json_result(json_result, args.output)

if __name__ == "__main__":
    main()
