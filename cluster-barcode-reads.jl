"""
cluster-barcode-reads.jl by Rohan Maddamsetti.

Let's cluster barcoding reads, allowing MAX_MIXMATCHES errors.
"""

using DataFrames, DataFramesMeta, CSV, BioSequences, BioAlignments

MAX_MISMATCHES = 10 #3

illumina_reads_df = CSV.read("ForwardRead_Cutoff3.csv", DataFrame)
readlength = unique([length(x) for x in illumina_reads_df.Sequence])
if length(readlength) != 1
    error("read length is not consistent: " * string(readlength))
else
    readlength = readlength[1]
end
similarity_threshold = round(1 - 3/readlength, sigdigits=2)

mutable struct Cluster
    ## the most common sequence is set as the reference.
    ref::LongDNASeq
    counts::Int64
end


# create an affine gap scoring model
affinegap = AffineGapScoreModel(
    match=1,
    mismatch=-1,
    gap_open=-1,
    gap_extend=-1
)

samples = unique(illumina_reads_df.Sample)
clusters_per_sample = []
for sample in samples
    sample_string = "Sample" * string(sample)    
    sample_df = @rsubset(illumina_reads_df, :Sample == sample)
    clusters = Vector{Cluster}()
    ## critical assumption: the sequences are read in sorted order,
    ## by the number of read counts.
    for (seq, count) in zip(sample_df.Sequence, sample_df.Count)
        DNAseq = LongDNASeq(seq)
        matched_cluster = false
        for cluster in clusters
            res = pairalign(GlobalAlignment(), DNAseq, cluster.ref, affinegap)
            aln = alignment(res)
            mismatches = count_mismatches(aln)
            insertions = count_insertions(aln)
            deletions = count_deletions(aln)
            if (mismatches + insertions + deletions) <= MAX_MISMATCHES
                matched_cluster = true
                cluster.counts += count
                break
            end
        end
        if matched_cluster == false
            newCluster = Cluster(DNAseq, count)
            push!(clusters, newCluster)
        end
    end
    push!(clusters_per_sample, clusters)
end

## write the clusters to file.
outfile = "merged_barcode_reads.csv"
outfh = open(outfile, "w")

header = "Sample,Sequence,Count"
println(outfh, header)
for (sample, clusters) in enumerate(clusters_per_sample)
    for cluster in clusters
        line = string(sample) * "," * string(cluster.ref) * "," * string(cluster.counts)
        println(outfh, line)
    end
end
close(outfh)
