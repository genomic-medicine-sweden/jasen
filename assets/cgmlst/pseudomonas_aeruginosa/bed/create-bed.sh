# Info about the reference genome can be found here:
# https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2
refseq_id=GCF_000006765.1
ncbi_ver=NC_002516.2
csvfile=${refseq_id}.csv
bedfile=${refseq_id}.bed
cgmlst_csv_url="https://www.cgmlst.org/ncs/schema/16115339/locus/?content-type=csv"

# Download CVS file with cgMLST schema for Pseudomonas aeruginosa
if [[ ! -f $csvfile ]]; then
    curl -o $csvfile $cgmlst_csv_url
    echo "Created: $csvfile"
else
    echo "Already exists: $csvfile"
fi;

# Convert the CSV table to BED format by:
# 1. Skipping first line
# 2. Shifting the start with -1 due to the 0-based BED format
# 3. Converting from sequence length to end-position
if [[ -f $bedfile ]]; then
    echo "Already exists: $bedfile"
else
    awk -F"\t" -v ncbi_ver=$ncbi_ver '(NR>1) { print ncbi_ver"\t" $4-1 "\t" $4-1+$5 }' $csvfile > $bedfile
    echo "Created: $bedfile"
fi;
