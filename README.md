# SlideVar
_Sliding Window Variation Identifier (SlideVar)_

SlideVar.pl is a computational pipeline for detecting lineage-specific changes based on sequence consistency. Its main method is to slide a window of size 1 over previously aligned sequences to assess sequence consistency at the single nucleotide level. It then searches for conserved positions that show lineage-specific variation in the target lineage. Because the pipeline only considers absolute character matches, it can also identify lineage-specific amino acid substitutions. Whether in genes or non-coding regulatory elements, functional sequences can be quite short (e.g. domains or transcription factor binding sites are often around 20bp). Therefore, SlideVar.pl is well suited for finding short segment specific changes within longer conserved regions. Overall, this pipeline uses windowed consistency scoring to identify lineage-specfic variations embedded in multiple sequence alignments.

---
## USAGE

This pipeline requires two input files: one is a multiple sequence alignment file in FASTA format, and the other is a marker text file used to label the target taxonomic group.


For a given multiple sequence alignment, SlideVar will by default select the first sequence as the reference sequence. For each of the other sequences, it will perform the single nucleotide resolution sequence consistency assessment within the given window size and output the consistency result for each species relative to the reference sequence.


Once all sequences have been processed, SlideVar will then perform lineage-specific mutation detection at each position of the sequences. If a position shows low conservation (e.g. match below 70%) in the marked species, but high similarity (e.g. match above 90%) in other species, it can be identified as a position with a lineage-specific change in the target group.



Usage:
    perl SlideVar.pl -in <input.fasta> -l <species.list> -w <window size> -con <conserved number> -div <changed number> -bn <background can not conserved species>

    -in input.fasta : sequence aligned in fasta format file
    -l species.list : species list with species marked
    -w window size : default 20 , the window size for calculating the conserved region
    -con conserved number : default 18 , identify as conserved region if the number of nucleotides in the window is large or equal to this value
    -div changed number : default 12 , identify as changed region if the number of nucleotides in the window is less or equal to this value
    -bn background can not conserved species : default 0 , how many species can be not conserved in the background species

    ---
    species list file format, one species per line and foreground species marked with '*':
    human
    mouse
    snake *
    frog
    caecilian *
    ...

    fasta file format: # species name should not contain '.' '-' '@' etc. ; '_' is allowed ; the first species should be the reference species
    >species_name
    AAGCTTGGG
    or
    >species_name.seqId
    AAGCTTGGG
