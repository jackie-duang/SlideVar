# SlideVar
_Sliding Window Variation Identifier (SlideVar)_

SlideVar.pl is a computational pipeline for detecting lineage-specific changes based on sequence consistency. Its main method is to slide a window of size 1 over previously aligned sequences to assess sequence consistency at the single nucleotide level. It then searches for conserved positions that show lineage-specific variation in the target lineage. Because the pipeline only considers absolute character matches, it can also identify lineage-specific amino acid substitutions. Whether in genes or non-coding regulatory elements, functional sequences can be quite short (e.g. domains or transcription factor binding sites are often around 20bp). Therefore, SlideVar.pl is well suited for finding short segment specific changes within longer conserved regions. Overall, this pipeline uses windowed consistency scoring to identify lineage-specfic variations embedded in multiple sequence alignments.

---
## USAGE


---
## Pipeline for finding lineage-specfic divergent Conserved Noncoding Elements


0. Make LAST database，for species with relatively distant phylogenetic relationships, MAM8 is recommended. (Here, we used LAST software to perform pairwised whole genome alignment. https://gitlab.com/mcfrith/last )
```
/path/to/last/bin/lastdb -u MAM8 ref ref.fa
```

1. Use lastal to get pairwised whole genome aligment, HOXD70 is better for distant phylogenetic relationships.
```
/path/to/last/bin/lastal -m10 -pHOXD70 /path/to/ref_MAM8_db query.fa > query.maf
```

2. Use last-split to get one-to-one alignment and sort the maf
```
/path/to/last/bin/last-split -m10 query.maf | /path/to/last/bin/last-split -r -m10 | /path/to/last/bin/maf-sort > query.sort.maf 
```

3. Combine all pairwised maf to multiple alignment maf
```
# remember to rename all pairwised alignment with maf format file to ref.query.sing.maf
/path/to/multiz/build/roast - T=. E=ref  "(ref, (species1, (species2 , species3)))" combined.maf > roast.maf.sh && sh roast.maf.sh
```
Or you can simply join all mafs with maf-join
```
/path/to/last/bin/maf-join *sort.maf > combined.maf # note this will only combine maf block which all species have sequence alignment
```

4. phastCons (optional) see: http://compgen.cshl.edu/phast/phastCons-HOWTO.html

5. You can also directly convert the maf block into fasta format and record tho position.
```
perl maf2fasta_by_speceis_list.pl <input.maf> <species.list> <min_length> <output_dir_prefix>
```

6. Run SlideVar.pl for each fasta file
```
for i in output_dir/*/*/Block*fasta; do 
	perl SlideVar.pl -in input.fasta -l species.list -w 20 -con 18 -div 12 -bn 2
done
```

7. Run mergeResults.pl to merge all results, the results in the file like "M1I2D3S4" means one site matched, 2 sites inserted, 3 sites deleted and 4 sites substituted. This also means the length of the window is 10 (1+2+3+4=10).
```
perl mergeResults.pl output_dir species.list output_file
```

---

```
perl SlideVar.pl -in <input.fasta> -l <species.list> -w <window size> -con <conserved number> -div <changed number> -bn <background diverged species>

    -in input.fasta : sequence aligned in fasta format file
    -l species.list : species list with species marked
    -w window size : default 20 , the window size for calculating the conserved region
    -con conserved number : default 18 , threshold for conserved region, identity >= conserved number is conserved
    -div diverged number : default 12 , threshold for divergent region, identity <= diverged number is divergent
    -bn  : default 0 , how many species could be divergent in background species (in case of some species are not conserved in background species because of assembly error or other reasons)

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

```
