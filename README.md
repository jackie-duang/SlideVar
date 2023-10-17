# SlideVar
Sliding Window Variation Identifier (SlideVar)

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
