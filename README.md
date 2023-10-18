# SlideVar
_Sliding Window Variation Identifier (SlideVar)_

SlideVar.pl is a computational pipeline for detecting lineage-specific changes based on sequence consistency. Its main method is to slide a window of size 1 over previously aligned sequences to assess sequence consistency at the single nucleotide level. It then searches for conserved positions that show lineage-specific variation in the target lineage. Because the pipeline only considers absolute character matches, it can also identify lineage-specific amino acid substitutions. Whether in genes or non-coding regulatory elements, functional sequences can be quite short (e.g. domains or transcription factor binding sites are often around 20bp). Therefore, SlideVar.pl is well suited for finding short segment specific changes within longer conserved regions. Overall, this pipeline uses windowed consistency scoring to identify lineage-specfic variations embedded in multiple sequence alignments.

---
## USAGE


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
