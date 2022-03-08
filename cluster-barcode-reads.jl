"""
cluster-barcode-reads.jl by Rohan Maddamsetti.

Cluster barcoding reads, allowing MAX_MIXMATCHES errors.
Also print out the location of matches found in (optional) reference FASTA sequences.

Input parameters:
=================
-- barcoding reads
-- maximum number of mismatches
-- reference FASTA sequences. One sequence per file, and multiple files are allowed.
-- path for output CSV file.

-- TODO: length of sequences to trim from 5' or 3' ends of the reads?


Assumed structure of barcode reads
==================================
Here is the structure of the sequencing data. The whole sequence is the primer, 
the first 24bp is the index, which should be removed. The next 74 bp is 5' sequence 
(it should be also constant among samples). Then the 29bp CRISPR repeat region is 
highly constant (gagttccccgcgccagcggggataaaccg), then followed with a 32bp spacer 
(the sequencing can only get the first 24bp, so all sequence after the CRISPR repeat 
can be used to find the spacer).

Output parameters:
==================
Count the occurrence of each barcode sequence in the reads, and find its 
matching location in the chromosome or plasmid.

Output is a CSV file with the following columns:

Sequence, Counts, Reference_Sequence, Reference_Location

The reference location is 1-indexed, and is an inclusive interval.

Usage examples:
==================
julia cluster-barcode-reads.jl ../data/first-test-data/ForwardRead_Cutoff3.csv

julia cluster-barcode-reads.jl ../data/second-test-data/test-samples.csv -r ../data/ref-seqs/mg1655-genome.fasta ../data/ref-seqs/a23-c-g8c1t-mcherry.fasta -o ../results/test_output.csv

"""

using DataFrames, DataFramesMeta, CSV, BioSequences, BioAlignments, FASTX, ArgParse


mutable struct Cluster
    ## the most common sequence is set as the reference.
    ref::LongDNASeq
    counts::Int64
end


function ClusterReads(illumina_reads_path, max_mismatches, affinegap)

    illumina_reads_df = CSV.read(illumina_reads_path, DataFrame)
    readlength = unique([length(x) for x in illumina_reads_df.Sequence])
    if length(readlength) != 1
        error("read length is not consistent: " * string(readlength))
    else
        readlength = readlength[1]
    end
    similarity_threshold = round(1 - 3/readlength, sigdigits=2)

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
                if (mismatches + insertions + deletions) <= max_mismatches
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

    return clusters_per_sample
end


function WriteClustersToFile(clusters_per_sample, outfile)
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
end


function get_CRISPR_spacer(readseq, CRISPR_repeat, affinegap, length_to_trim = 95)
    """
    Assumed structure of barcode reads
    ==================================
    Here is the structure of the sequencing data. The whole sequence is the primer, 
    the first 24bp is the index, which should be removed. The next 74 bp is 5' sequence 
    (it should be also constant among samples). Then the 29bp CRISPR repeat region is 
    highly constant (gagttccccgcgccagcggggataaaccg), then followed with a 32bp spacer 
    (the sequencing can only get the first 24bp, so all sequence after the CRISPR repeat 
    can be used to find the spacer).
    """
    trimmed_read = readseq[length_to_trim:end]
    ## align the trimmed read to the CRISPR repeat.
    res = pairalign(GlobalAlignment(), trimmed_read, CRISPR_repeat, affinegap)
    aln = alignment(res)    
    ## get the last position of the repeat; the spacer follows.
    spacer_start = ref2seq(aln, lastindex(CRISPR_repeat))[1] + 1
    spacer = trimmed_read[spacer_start:end]
    return spacer
end


function GetPairwiseAlignmentSeqStart(aln)
    seq = aln.a.seq
    anchors = aln.a.aln.anchors

    for 
    
    seqpos = anchors[1].seqpos
    
seq = aln.a.seq
    ref = aln.b
    anchors = aln.a.aln.anchors
    # width of position numbers
    posw = ndigits(max(anchors[end].seqpos, anchors[end].refpos)) + 1

    i = 0

    refpos = anchors[1].refpos

function AlignClustersToRefSeq(clusters_per_sample, refseq, affinegap, CRISPR_repeat)

    for (sample, clusters) in enumerate(clusters_per_sample)
        for cluster in clusters
            matched_reference = false

            spacer = get_CRISPR_spacer(cluster.ref, CRISPR_repeat, affinegap)
            res = pairalign(LocalAlignment(), spacer, refseq, affinegap)
            aln = alignment(res)
            println(aln)


            println(first(aln.anchors).seqpos)
            last(aln.anchors).seqpos)
            first(aln.anchors).refpos, '-', last(aln.anchors).refpos)
            
        end
    end
    
end


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "barcode_read_csv_path"
        help = "csv file of barcoded reads and raw counts."
        arg_type = String
        required = true
        "max_mismatches"
        help = "maximum number of mismatches allowed during clustering"
        arg_type = Int64
        default = 10
        "--reference_fasta_paths", "-r"
        help = "reference fasta files, for finding CRISPR spacer source."
        nargs = '+'
        arg_type = String
        "--outfile", "-o"
        help = "outfile for the clustered reads. default is merged_barcode_reads.csv."
        arg_type = String
        default = "../results/merged_barcode_reads.csv"
    end
    return parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    
    illumina_reads_path = parsed_args["barcode_read_csv_path"]
    max_mismatches = parsed_args["max_mismatches"]
    reference_fasta_paths = parsed_args["reference_fasta_paths"]
    outfile = parsed_args["outfile"]

    ## create an affine gap scoring model
    affinegap = AffineGapScoreModel(
        match=1,
        mismatch=-1,
        gap_open=-1,
        gap_extend=-1
    )
    
    ## This CRISPR repeat is constant, defined by Yi.
    CRISPR_repeat = dna"gagttccccgcgccagcggggataaaccg"

    ## cluster the reads.
    clusters_per_sample = ClusterReads(illumina_reads_path, max_mismatches, affinegap)
    ## align the clusters to each reference FASTA sequence, if they exist.

    for ref_fasta_file in reference_fasta_paths
        reader = open(FASTA.Reader, ref_fasta_file)
        for record in reader
            refseq = sequence(record)
            ## now, align each read to the current reference.
            AlignClustersToRefSeq(clusters_per_sample, refseq, affinegap, CRISPR_repeat)
        end
        close(reader)
    end
    
    ## write the clusters to file.
    ##WriteClustersToFile(clusters_per_sample, outfile)

end


main()
