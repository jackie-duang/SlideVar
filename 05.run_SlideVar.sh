for i in output_dir/*/*/Block*fasta; do 
	perl SlideVar.pl -in input.fasta -l species.list -w 20 -con 18 -div 12 -bn 2
done
