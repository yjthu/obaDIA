#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

sub usage
{
        print STDERR <<USAGE;
=========================================================================
Description		protein functional enrichment pipeline

Options
		-directory-kobas-output		<s> : output files directory
		-samplename 			<s> : sample name for output files
		-species 			<s> : species for kegg or background annotate file
		-diffGeneSeq 			<s> : protein sequence fasta
		-goanno 			<s> : goanno file
		-diffid 			<s> : protein id list
		-updown 			<s> : updown file
		-difftb 			<s> : diffrential expressed protein information file
		-ko 				<s> : if or not ko enrichment [yes/no]
		-spe_go 			<s> : species for kobas GO and Rectome [na/spe]
		-pfanno 			<s> : pfam_anno file
		-summary 			<s> : protein length file
		-bg 				<s> : background file prefix[*.spe.annotate]
		-mod				<s> : mode for backgroud, T for -spe and E for -bg(default)
		-h|?|help   			    :  Show this help
=========================================================================
USAGE
}

my $base_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA';
my $db_dir='/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db';
my $bin = "$base_dir/bin/Enrich";
my $gotype = "$db_dir/go.type";
my $kobasDB = "/storage/data/PUBLIC/databases/KOBAS_3.0_db";

my $Rscript = "Rscript";
my $diamond="diamond";

my ($outdir,$samplename,$species,$diffseq,$goanno,$diffid,$updown,$difftb,$ko,$spe_go,$pfanno,$summary,$bg,$mod,$help);
GetOptions (
	"directory-kobas-output=s" => \$outdir,
	"samplename=s"=>\$samplename,
	"species=s" => \$species,
	"diffGeneSeq=s" => \$diffseq,
	"goanno=s" => \$goanno,
	"diffid=s" => \$diffid,
	"updown=s" => \$updown,
	"difftb=s" => \$difftb,
	"ko=s" => \$ko,
	"spe_go=s" => \$spe_go,
	"pfanno=s" => \$pfanno,
	"summary=s" => \$summary,
	"bg=s" => \$bg,
	"mod=s" => \$mod,
	"h|?|help"=>\$help,
);


if(!defined($outdir) || !defined($samplename) || !defined($species) || !defined($diffseq) || defined($help)){
        &usage;
        exit 0;
}



unless (-d $outdir){
	!system "mkdir -p $outdir" or die "Something went wrong when mkdir for $outdir!";
}

if(-e "$kobasDB/seq_pep/$species.pep.fasta.gz" || -e "$kobasDB/sqlite3/$species.db.gz"){
	`gunzip $kobasDB/seq_pep/$species.pep.fasta.gz`;
	`gunzip $kobasDB/sqlite3/$species.db.gz`;
}
if(!-e  "$kobasDB/seq_pep/$species.dmnd"){
	`$diamond makedb --in $kobasDB/seq_pep/$species.pep.fasta -d $kobasDB/seq_pep/$species`;
	`chmod 777 $kobasDB/seq_pep/$species.dmnd`;
}

if($spe_go ne $species && $spe_go ne 'na'){
	if(-e "$kobasDB/seq_pep/$spe_go.pep.fasta.gz" || -e "$kobasDB/sqlite3/$spe_go.db.gz"){
		`gunzip $kobasDB/seq_pep/$spe_go.pep.fasta.gz`;
		`gunzip $kobasDB/sqlite3/$spe_go.db.gz`;
	}
	if(!-e  "$kobasDB/seq_pep/$spe_go.dmnd"){
		`$diamond makedb --in $kobasDB/seq_pep/$spe_go.pep.fasta -d $kobasDB/seq_pep/$spe_go`;
		`chmod 777 $kobasDB/seq_pep/$spe_go.dmnd`;
	}
}

$ko ||= 'no';
$spe_go ||= $species;
$mod ||= 'E';
my $identify_bg = $species;

my $cmd='';

# kobas kegg enrich
if($species eq 'ko'){
	$cmd .= "\n\n$diamond blastp -d $kobasDB/seq_pep/$species.dmnd -q $diffseq -o $outdir/$samplename.$species.dmout --outfmt 6 --tmpdir . -k 1 --quiet\n";
	$cmd .= "kobas-annotate -i $outdir/$samplename.$species.dmout -t blastout:tab -s $species -o $outdir/$samplename.$species.annotate\n";
	$cmd .= "kobas-identify -f $outdir/$samplename.$species.annotate -b $outdir/$samplename.$species.annotate -d K -o $outdir/$samplename.identify.Kegg\n";
	$cmd .= "perl $bin/add_id.pl $outdir/$samplename.$species.annotate $outdir/$samplename.identify.Kegg >$outdir/$samplename.KeggEnrich\n";
	$spe_go = 'na';
}

if($species ne 'ko'){
	if($bg){
		if($mod eq 'E'){
			$identify_bg="$bg.$species.annotate";
		}
	}
	$cmd .= "$diamond blastp -d $kobasDB/seq_pep/$species.dmnd -q $diffseq -o $outdir/$samplename.$species.dmout --outfmt 6 --tmpdir . --more-sensitive -k 1 --quiet\n";
	$cmd .= "kobas-annotate -i $outdir/$samplename.$species.dmout -t blastout:tab -s $species -o $outdir/$samplename.$species.annotate\n";
	$cmd .= "kobas-identify -f $outdir/$samplename.$species.annotate -b $identify_bg -d K -o $outdir/$samplename.identify.Kegg\n";
	$cmd .= "perl $bin/add_id.pl $outdir/$samplename.$species.annotate $outdir/$samplename.identify.Kegg >$outdir/$samplename.KeggEnrich\n";
	$cmd .= "perl $bin/extract_top20.pl $outdir/$samplename.KeggEnrich >$outdir/$samplename.KeggEnrich.top20\n";
	$cmd .= "$Rscript $bin/Pathwayscatter.R $samplename  $outdir $outdir/$samplename.KeggEnrich.top20\n\n";
}

if($spe_go ne 'na'){
	if($spe_go ne $species){
		$cmd .= "$diamond blastp -d $kobasDB/seq_pep/$spe_go.dmnd -q $diffseq -o $outdir/$samplename.$spe_go.dmout --outfmt 6 --tmpdir . --more-sensitive -k 1 --quiet\n";
		$cmd .= "kobas-annotate -i $outdir/$samplename.$spe_go.dmout -t blastout:tab -s $spe_go -o $outdir/$samplename.$spe_go.annotate\n";
	}
	
	# kobas go enrich
	$identify_bg = $spe_go;
	if($bg){
		if($mod eq 'E'){
			$identify_bg="$bg.$spe_go.annotate";
		}
	}
	$cmd .= "kobas-identify -f $outdir/$samplename.$spe_go.annotate -b $identify_bg -d G -o $outdir/$samplename.identify.GO\n\n";
	if($goanno && $gotype){
		$cmd .= "perl $bin/go_resource.pl $goanno $outdir/go\n";
		$cmd .= "grep -v '#' $outdir/$samplename.identify.GO|grep -v '^-'|grep -v '^\$'|perl $bin/add_go_type.pl $gotype - >$outdir/$samplename.identify.GO.xls\n";
		$cmd .= "$Rscript $bin/goBar.R $outdir/$samplename.identify.GO.xls $outdir $samplename\n";
		if($diffid){
			$cmd .= "$Rscript $bin/topGO.R $outdir $outdir/go.resource $outdir/$samplename.identify.GO.xls $diffid $outdir/$samplename\n\n";
		}
	}
	
	# kobas rectome enrich
	$cmd .= "kobas-identify -f $outdir/$samplename.$spe_go.annotate -b $identify_bg -d R -o $outdir/$samplename.identify.Rectome\n\n";
	$cmd .= "grep -v '#' $outdir/$samplename.identify.Rectome|grep -v '^-'|grep -v '^\$'|sed '1i#Term\\tDatabase\\tID\\tInput_number\\tBackground_number\\tP_Value\\tCorrected_P_Value\\tInput\\tHyperlink' >$outdir/$samplename.identify.Rectome.xls\n";
	$cmd .= "perl $bin/extract_top20.pl $outdir/$samplename.identify.Rectome >$outdir/$samplename.identify.Rectome.top20\n";
	$cmd .= "Rscript $bin/Pathwayscatter2.R $samplename.Rectome $outdir $outdir/$samplename.identify.Rectome.top20\n\n";
}


# goseq pfam enrich
if($pfanno && $summary && $diffid){
	$cmd .= "$bin/runPfamEnrich.sh $summary $pfanno $diffid $samplename $outdir\n\n";
}

# goseq go enrich
if($spe_go eq 'na'){
	if($goanno && $summary && $diffid){
		$cmd .= "$bin/runGOEnrich.sh $summary $goanno $diffid $samplename $outdir\n\n";
	}
}

# kegg web 
if($updown){
	if($bg && $difftb){
		$cmd .= "perl $bin/kobasMapColor.v3.pl $bg.$species.annotate $updown $bg.KeggEnrich $outdir/$samplename.KeggEnrich $outdir/$samplename.annomap >$outdir/$samplename.KeggEnrich.xls\n";
		$cmd .= "python $bin/kegg_pathway_web3.py --table $outdir/$samplename.KeggEnrich.xls --diff $updown --diff2 $difftb --anno $outdir/$samplename.annomap\n\n";
	}else{
		$cmd .= "perl $bin/kobasMapColor.v2.pl $outdir/$samplename.$species.annotate $updown $outdir/$samplename.KeggEnrich >$outdir/$samplename.KeggEnrich.xls\n";
		$cmd .= "python $bin/kegg_pathway_web2.py --table $outdir/$samplename.KeggEnrich.xls --diff $updown\n\n";
	}
}else{
	$cmd .= "perl $bin/kobasMapColor.pl $outdir/$samplename.$species.annotate $outdir/$samplename.KeggEnrich >$outdir/$samplename.KeggEnrich.xls\n";
	$cmd .= "python $bin/kegg_pathway_web.py --table $outdir/$samplename.KeggEnrich.xls\n\n";
}

if($ko eq 'yes'){
	$cmd .= "kobas-identify -f $outdir/$samplename.$species.annotate -b $bg.$species.annotate -d K -o $outdir/$samplename.identify.Kegg\n";
	$cmd .= "perl $bin/add_id.pl $outdir/$samplename.$species.annotate $outdir/$samplename.identify.Kegg >$outdir/$samplename.KeggEnrich\n";
	$cmd .= "perl $bin/extract_top20.pl $outdir/$samplename.KeggEnrich >$outdir/$samplename.KeggEnrich.top20\n";
	$cmd .= "$Rscript $bin/Pathwayscatter.R $samplename $outdir $outdir/$samplename.KeggEnrich.top20\n\n";
	if($updown){
		if($bg && $difftb){
			$cmd .= "perl $bin/kobasMapColor.v3.pl $bg.$species.annotate $updown $bg.KeggEnrich $outdir/$samplename.KeggEnrich $outdir/$samplename.annomap >$outdir/$samplename.KeggEnrich.xls\n";
			$cmd .= "python $bin/kegg_pathway_web3.py --table $outdir/$samplename.KeggEnrich.xls --diff $updown --diff2 $difftb --anno $outdir/$samplename.annomap\n\n";
		}else{
			$cmd .= "perl $bin/kobasMapColor.v2.pl $outdir/$samplename.$species.annotate $updown $outdir/$samplename.KeggEnrich >$outdir/$samplename.KeggEnrich.xls\n";
			$cmd .= "python $bin/kegg_pathway_web2.py --table $outdir/$samplename.KeggEnrich.xls --diff $updown\n\n";
		}
	}else{
		$cmd .= "perl $bin/kobasMapColor.pl $outdir/$samplename.$species.annotate $outdir/$samplename.KeggEnrich >$outdir/$samplename.KeggEnrich.xls\n";
		$cmd .= "python $bin/kegg_pathway_web.py --table $outdir/$samplename.KeggEnrich.xls \n\n";
	}
}

open OUT, ">$outdir/.$samplename.enrichP.sh" or die $!;
print OUT "$cmd";
close OUT;
`chmod +x $outdir/.$samplename.enrichP.sh`;
