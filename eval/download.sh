if [ ! -r ecoli.fa ]
then
  wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz
  gunzip Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz
  mv Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa ecoli.fa
fi

if [ ! -r celegans.fa ]
then
  wget ftp://ftp.ensembl.org/pub/release-98/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
  gunzip Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz
  mv Caenorhabditis_elegans.WBcel235.dna.toplevel.fa celegans.fa
fi

WORKINGDIR=`pwd`

if [ $WORKINGDIR == '/work-zfs/mschatz1/mkirsche/sapling/genomes' ] && [ ! -r human.fa ]
then
  wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
  gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
  mv GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta human.fa
fi

if [ $WORKINGDIR == '/work-zfs/mschatz1/mkirsche/sapling/genomes' ] && [ ! -r wheat.fa ]
then
  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-45/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
  gunzip Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
  mv Triticum_aestivum.IWGSC.dna.toplevel.fa wheat.fa
fi

