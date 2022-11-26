"""
cluster-csv-barcode-reads.jl by Rohan Maddamsetti.

This code assumes that the input reads are in CSV format.

Cluster barcoding reads, allowing MAX_MIXMATCHES errors.
Also print out the location of matches found in (optional) reference FASTA sequences.

Input parameters:
=================
-- barcoding reads in CSV format.
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


function get_CRISPR_spacer(readseq, CRISPR_repeat, affinegap, length_to_trim = 95, max_CRISPR_mismatches = 2)
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
    ## error checking: if the best alignment poorly matches the CRISPR repeat,
    num_CRISPR_matches = count_matches(aln)
    if (num_CRISPR_matches < (length(CRISPR_repeat) - max_CRISPR_mismatches))
        spacer = "" ## then return an empty spacer string.
    else
        ## get the last position of the repeat; the spacer follows.
        spacer_start = ref2seq(aln, lastindex(CRISPR_repeat))[1] + 1
        spacer = trimmed_read[spacer_start:end]
    end
    return spacer
end


function GetPairwiseAlignmentSeqStart(pairaln::PairwiseAlignment)
    return pairaln.a.aln.anchors[1].seqpos + 1
end


function GetPairwiseAlignmentRefStart(pairaln::PairwiseAlignment)
    return pairaln.a.aln.anchors[1].refpos + 1
end


function GetPairwiseAlignmentSeqEnd(pairaln::PairwiseAlignment)
    return pairaln.a.aln.anchors[end].seqpos
end


function GetPairwiseAlignmentRefEnd(pairaln::PairwiseAlignment)
    return pairaln.a.aln.anchors[end].refpos
end


function AlignClusterToRefSeq(cluster, refseq, affinegap, CRISPR_repeat)
    spacer = get_CRISPR_spacer(cluster.ref, CRISPR_repeat, affinegap)
    ## error checking: if the CRISPR repeat does not exist in the
    ## reference sequence, then an empty spacer is returned.
    if (length(spacer) == 0) ## in the case of junk reads, return nonsense values.
        refstart = -1
        refend = -1
        aln_score = -1
    else
        res = pairalign(LocalAlignment(), spacer, refseq, affinegap)
        aln_score = score(res)
        aln = alignment(res)            
        refstart = GetPairwiseAlignmentRefStart(aln)
        refend = GetPairwiseAlignmentRefEnd(aln)
    end
    return (spacer, refstart, refend, aln_score)
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


function WriteAlignedClustersToFile(clusters_per_sample, reference_fasta_paths, affinegap, CRISPR_repeat, outfile, reference_alignment_threshold = 20)
    outfh = open(outfile, "w")
    header = "Sample,Sequence,Count,Spacer,Reference,ReferenceStart,ReferenceEnd,AlignmentScore"
    println(outfh, header)

    ## align to each record in each fasta file. 
    for ref_fasta_file in reference_fasta_paths
        reader = open(FASTA.Reader, ref_fasta_file)
        for record in reader
            record_name = FASTA.identifier(record) * "_" * FASTA.description(record)
            refseq = sequence(record)
            ## now, align each read to the current reference.
            for (sample, clusters) in enumerate(clusters_per_sample)
                for cluster in clusters
                    spacer, refstart, refend, aln_score = AlignClusterToRefSeq(cluster, refseq, affinegap, CRISPR_repeat)
                    ## error handling.
                    ## case 1: CRISPR repeat is not in the read.
                    if (length(spacer) == 0)
                        continue
                    end
                    ## case 2: score threshold for the alignment to reference is too low.
                    if (aln_score < reference_alignment_threshold)
                        continue
                    end
                    
                    ## write the alignment data to file.
                    fields = [string(sample), string(cluster.ref), string(cluster.counts), string(spacer), record_name, refstart, refend, aln_score]
                    line = join(fields, ",")
                    println(outfh, line)
                end
            end
        end
        ## close the current reference FASTA file handle.
        close(reader)
    end
    ## close the output file handle.
    close(outfh)
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
        default = 5
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

    ## if reference fasta files have been supplied, then align and write to file.
    if (length(reference_fasta_paths) > 0)
        WriteAlignedClustersToFile(clusters_per_sample, reference_fasta_paths, affinegap, CRISPR_repeat, outfile)
    else
        ## otherwise, just write the clusters to file.
        WriteClustersToFile(clusters_per_sample, outfile)
    end
    ## end of the main function.
end


main()
