"""
cluster-fastq-barcode-reads.jl by Rohan Maddamsetti.

This code assumes that the input reads are in FASTQ.GZ format
(compressed FASTQ reads).

This code also assume that ONLY forward reads from paired-end sequencing are provided.

Cluster barcoding reads, allowing MAX_MIXMATCHES errors.
Also print out the location of matches found in (optional) reference FASTA sequences.

Input parameters:
=================
-- barcoding reads in *.fastq.gz format
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
julia cluster-fastq-barcode-reads.jl ../data/third-test-data/

julia cluster-barcode-reads.jl ../data/third-test-data/ -r ../data/ref-seqs/a27-c-orit-g8c1t-mcherry.fasta -o ../results/genewiz_forward_merged_barcodes.csv

julia cluster-barcode-reads.jl ../data/third-test-data/ -r ../data/ref-seqs/a27-c-orit-g8c1t-mcherry.fasta -o ../results/genewiz_forward_merged_barcodes.csv

julia cluster-barcode-reads.jl ../data/2022-11-16-merged-data/ -r ../data/ref-seqs/a27-c-orit-g8c1t-mcherry.fasta ../data/ref-seqs/c16-c-a-train-mcherry.fasta -o ../results/merged_barcodes_2022-11-16.csv

"""

using DataFrames, DataFramesMeta, CSV, BioSequences, BioAlignments, FASTX, CodecZlib, ArgParse


mutable struct Cluster
    ## the most common sequence is set as the reference.
    ref::LongDNASeq
    counts::Int64
end


function TrimAndCountReads(fastq_gz_reads_dir, length_to_trim = 188, length_to_keep = 32, full_length=250)
    
    files_in_reads_dir = readdir(fastq_gz_reads_dir, join=true)
    fastq_gz_files = [x for x in files_in_reads_dir if endswith(x,".fastq.gz")]

    sample_to_sequence_count_dict = Dict()
    ## update the dictionary using the fastq_files in the directory.
    for fastq_gz_f in fastq_gz_files
        ## the name of the sample comes from the fastq.gz filename.
        sample = string(split(split(basename(fastq_gz_f),".fastq.gz")[1], "_")[1])

        if !haskey(sample_to_sequence_count_dict, sample)
            sample_to_sequence_count_dict[sample] = Dict()
        end

        for record in FASTQ.Reader(GzipDecompressorStream(open(fastq_gz_f)))
            my_seq = sequence(record)
            ## skip reads that are not full_length.
            if length(my_seq) < full_length
                continue
            end
            my_trimmed_seq = my_seq[length_to_trim:length_to_trim+length_to_keep]
            if !haskey(sample_to_sequence_count_dict[sample], my_trimmed_seq)
                sample_to_sequence_count_dict[sample][my_trimmed_seq] = 1
            else
                sample_to_sequence_count_dict[sample][my_trimmed_seq] += 1
            end
        end        
    end

    illumina_reads_df = DataFrame(Sample=String[], Sequence=String[], Count=Int64[])
    ## now turn the dictionary into a dataframe.
    for (sample, seq_count_dict) in sample_to_sequence_count_dict
        for (seq, seq_count) in seq_count_dict
            push!(illumina_reads_df, [sample seq seq_count])
        end
    end
    
    ## CRITICAL: sort the sequences by sample and the number of read counts.
    sort!(illumina_reads_df, [:Sample, :Count], rev=true)

    return illumina_reads_df
end


##function TrimAndCountForwardReads(fastq_gz_reads_dir, length_to_trim = 188, length_to_keep = 32, full_length=250)
##    return TrimAndCountReads(fastq_gz_reads_dir, length_to_trim = 188, length_to_keep = 32, full_length=250)


##function TrimAndCountReverseReads(fastq_gz_reads_dir, length_to_trim = 188, length_to_keep = 32, full_length=250)
##    return TrimAndCountReads(fastq_gz_reads_dir, length_to_trim = 188, length_to_keep = 32, full_length=250)



function ClusterReads(illumina_reads_df, max_mismatches, affinegap)
    
    readlength = unique([length(x) for x in illumina_reads_df.Sequence])
    if length(readlength) != 1
        error("read length is not consistent: " * string(readlength))
    else
        readlength = readlength[1]
    end
    similarity_threshold = round(1 - 3/readlength, sigdigits=2)

    samples = unique(illumina_reads_df.Sample)
    clusters_per_sample = Dict()
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
        clusters_per_sample[sample] = clusters
    end

    return clusters_per_sample
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


function AlignClusterToRefSeq(cluster, refseq, affinegap)    
    res = pairalign(LocalAlignment(), cluster.ref, refseq, affinegap)
    aln_score = score(res)
    aln = alignment(res)            
    refstart = GetPairwiseAlignmentRefStart(aln)
    refend = GetPairwiseAlignmentRefEnd(aln)
    return (refstart, refend, aln_score)
end


function WriteClustersToFile(clusters_per_sample, outfile)
    outfh = open(outfile, "w")
    header = "Sample,Sequence,Count"
    println(outfh, header)
    for (sample, clusters) in clusters_per_sample
        for cluster in clusters
            line = string(sample) * "," * string(cluster.ref) * "," * string(cluster.counts)
            println(outfh, line)
        end
    end
    close(outfh)
end


function WriteAlignedClustersToFile(clusters_per_sample, reference_fasta_paths, affinegap, outfile, reference_alignment_threshold = 20)
    outfh = open(outfile, "w")
    header = "Sample,Sequence,Count,Reference,ReferenceStart,ReferenceEnd,AlignmentScore"
    println(outfh, header)

    ## align to each record in each fasta file. 
    for ref_fasta_file in reference_fasta_paths
        reader = open(FASTA.Reader, ref_fasta_file)
        for record in reader
            record_name = FASTA.identifier(record) * "_" * FASTA.description(record)
            refseq = sequence(record)
            ## now, align each read to the current reference.
            for (sample, clusters) in clusters_per_sample
                for cluster in clusters
                    refstart, refend, aln_score = AlignClusterToRefSeq(cluster, refseq, affinegap)

                    ## skip if score threshold for the alignment to reference is too low.
                    if (aln_score < reference_alignment_threshold)
                        continue
                    end
                    
                    ## write the alignment data to file.
                    fields = [string(sample), string(cluster.ref), string(cluster.counts), record_name, refstart, refend, aln_score]
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
        "barcode_read_fastq_gz_dir"
        help = "directory containing barcoded reads in fastq.gz format."
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
    
    fastq_gz_reads_dir = parsed_args["barcode_read_fastq_gz_dir"]
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
    
    ## trim and count the reads.
    illumina_reads_df = TrimAndCountReads(fastq_gz_reads_dir)
    
    ## cluster the reads.
    ## This function takes about 45 seconds for one sample.
    clusters_per_sample = ClusterReads(illumina_reads_df, max_mismatches, affinegap)
    
    ## if reference fasta files have been supplied, then align and write to file.
    if (length(reference_fasta_paths) > 0)
        WriteAlignedClustersToFile(clusters_per_sample, reference_fasta_paths, affinegap, outfile)
    else
        ## otherwise, just write the clusters to file.
        WriteClustersToFile(clusters_per_sample, outfile)
    end
    ## end of the main function.
end


main()
