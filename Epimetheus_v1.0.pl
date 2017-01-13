=head
Epimetheus Version 1.0
Last modified March 27, 2015

For bugs and suggestions, please contact mohameds@igbmc.fr
=cut


use IO::File;
use File::Basename;
use Getopt::Long;

$tooldir = dirname (__FILE__);

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}


sub wig2bed {
foreach $hist (sort keys %fileList)
{
	$time=localtime();
	print STDOUT "$hist: Generating BED file from normalized WIG file with Read Lenght:$rl Fragment Length:$ext [$time]\n";
	foreach $samp (sort keys %{$fileList{$hist}})
	{
		$fileCur=$fileList{$hist}{$samp};
		$id="$hist-$samp";$base="$opFolder/$hist/$samp/$id";
		$opZ=`file $fileCur`;
		if($opZ=~/compressed data/)
		{
			$cmd="zcat $fileCur >$base.raw.bed";
			`$cmd`;
			$fileCur="$base.raw.bed";
		}
		$rl=$readLenOut;
		$ext=$extReads;
		print STDOUT "Processing $samp\tCurrent File:$fileCur Normalized BED:$base.norm.bed\n";
		open(BED,$fileCur);
		while(<BED>)
		{
			chomp();
			@a=split("\t",$_);
			if($a[5]=~/\+/){$mid=$a[1]+($ext/2);}else{$mid=$a[2]-($ext/2);}
			$l1=$mid+$bin-($mid) % $bin if($mid % $bin);
			if($a[1] % $bin){
				$l1=($a[1]+$bin-($a[1]) % $bin)-$bin;
			}
			else{
			$l1=$a[1];
			}
			if($keepClonalReads==0){$bedF{$a[0]}{$l1}{$a[1]}{$a[5]}=1;}
			else{$bedF{$a[0]}{$l1}{$a[1]}{$a[5]}++;}
			$bedC{$a[0]}{$l1}++;
		}
		close BED;
		open(RAW,"$base.raw.bedgraph") || "Couldn't open file $base.raw.bedgraph\n";
		open(NORM,"$base.norm.bedgraph") || "Couldn't open file $base.norm.bedgraph\n";
		open(OUT,">$base.norm.bed") || "Couldn't open file $base.norm.bed\n";
		$chr="chr0";
		while(<NORM>)
		{
			chomp();
			$rawB=0;
			chomp($rin=<RAW>);
			@a=split("\t",$_);
			$bin1=int($a[3]+0.5);
			@ra=split("\t",$rin);
			$rawB=int($ra[3]+0.5);
			$mkk++;
			if($keepClonalReads==1){
				if($bin1==$rawB)
				{
					foreach $k (sort keys %{$bedF{$a[0]}{$a[1]}})
					{
						foreach $m (sort keys %{$bedF{$a[0]}{$a[1]}{$k}})
						{
							for($i=0;$i<$bedF{$a[0]}{$a[1]}{$k}{$m};$i++)
							{
								if($m=~/\+/){print OUT "$a[0]\t",$k-($ext/2),"\t",$k-($ext/2)+$rl,"\tseq1\t0\t$m\n";}else{print OUT "$a[0]\t",$k+($ext/2)-$rl,"\t",$k+($ext/2),"\tseq1\t0\t$m\n";}
							}
						}
					}
				}
				elsif($bin1<$rawB)
				{
					$diff=$bin1;
					foreach $k (keys %{$bedF{$a[0]}{$a[1]}})
					{
						foreach $m (keys %{$bedF{$a[0]}{$a[1]}{$k}})
						{
							for($i=0;$i<$bedF{$a[0]}{$a[1]}{$k}{$m};$i++)
							{
								$diff--;
								if($diff==0){last;}
								if($m=~/\+/){print OUT "$a[0]\t",$k-($ext/2),"\t",$k-($ext/2)+$rl,"\tseq2\t0\t$m\n";}else{print OUT "$a[0]\t",$k+($ext/2)-$rl,"\t",$k+($ext/2),"\tseq2\t0\t$m\n";}
							}
						if($diff==0){last;}
						}
					if($diff==0){last;}
					}
				}
				elsif($rawB<$bin1)
				{
					$diff=$bin1-$rawB;
					foreach $k (sort keys %{$bedF{$a[0]}{$a[1]}})
					{
						foreach $m (sort keys %{$bedF{$a[0]}{$a[1]}{$k}})
						{
							for($i=0;$i<$bedF{$a[0]}{$a[1]}{$k}{$m};$i++)
							{
								if($m=~/\+/){print OUT "$a[0]\t",$k-($ext/2),"\t",$k-($ext/2)+$rl,"\tseq3\t0\t$m\n";}else{print OUT "$a[0]\t",$k+($ext/2)-$rl,"\t",$k+($ext/2),"\tseq3\t0\t$m\n";}
							}
						}
					}
					$p=0;$x=-1;
					for($i=0;$i<$diff;)
					{
						$p=0;
						while($p<$bin){$x = $a[1] + int(rand($a[2] - $a[1]));if(exists $bedF{$a[0]}{$a[1]}{$x}{'+'} && exists $bedF{$a[0]}{$a[1]}{$x}{'-'}){$p++;$x=-1;}else{last;}}
						if($x<0){
							while($i<$diff){$x = $a[1] + int(rand($a[2] - $a[1]));$st=int(rand()+0.5);if($st==1){print OUT "$a[0]\t",$x-($ext/2),"\t",$x-($ext/2)+$rl,"\tseq3\t0\t+\n";}else{print OUT "$a[0]\t",$x+($ext/2)-$rl,"\t",$x+($ext/2),"\tseq3\t0\t+\n";}$i++;}
						}
						$st=int(rand()+0.5);
						if(exists $bedF{$a[0]}{$a[1]}{$x}{$st}){if($st==1){$st=0;}else{$st=1;}}
						if($st==1){
						$bedF{$a[0]}{$a[1]}{$x}{$st}++;
						print OUT "$a[0]\t",$x-($ext/2),"\t",$x-($ext/2)+$rl,"\tseq3\t0\t+\n";
						$i++;
						}
						else{
						$bedF{$a[0]}{$a[1]}{$x}{$st}++;
						print OUT "$a[0]\t",$x-($ext/2),"\t",$x-($ext/2)+$rl,"\tseq3\t0\t+\n";
						$i++;
						}
					}
				}
			}
			else
			{
				if($bin1==$rawB)
				{
					foreach $k (sort keys %{$bedF{$a[0]}{$a[1]}})
					{
						foreach $kk (sort keys %{$bedF{$a[0]}{$a[1]}{$k}})
						{
							if($kk=~/\+/){print OUT "$a[0]\t",$k-($ext/2),"\t",$k-($ext/2)+$rl,"\tseq1\t0\t$kk\n";}else{print OUT "$a[0]\t",$k+($ext/2)-$rl,"\t",$k+($ext/2),"\tseq1\t0\t$kk\n";}
						}
					}
				}
				elsif($bin1<$rawB)
				{
					$diff=$bin1;
					foreach $k (keys %{$bedF{$a[0]}{$a[1]}})
					{
						foreach $kk (sort keys %{$bedF{$a[0]}{$a[1]}{$k}})
						{
							$diff--;
							if($diff==0){last;}
							if($kk=~/\+/){print OUT "$a[0]\t",$k-($ext/2),"\t",$k-($ext/2)+$rl,"\tseq1\t0\t$kk\n";}else{print OUT "$a[0]\t",$k+($ext/2)-$rl,"\t",$k+($ext/2),"\tseq1\t0\t$kk\n";}
						}
						if($diff==0){last;}
					}
				}
				elsif($rawB<$bin1)
				{
					$diff=$bin1-$rawB;
					foreach $k (sort keys %{$bedF{$a[0]}{$a[1]}})
					{
						foreach $m (sort keys %{$bedF{$a[0]}{$a[1]}{$k}})
						{
							if($m=~/\+/){print OUT "$a[0]\t",$k-($ext/2),"\t",$k-($ext/2)+$rl,"\tseq3\t0\t$m\n";}else{print OUT "$a[0]\t",$k+($ext/2)-$rl,"\t",$k+($ext/2),"\tseq3\t0\t$m\n";}
						}
					}
					$p=0;$x=-1;
					for($i=0;$i<$diff;)
					{
						$p=0;
						while($p<$bin){$x = $a[1] + int(rand($a[2] - $a[1]));if(exists $bedF{$a[0]}{$a[1]}{$x}{'+'} && exists $bedF{$a[0]}{$a[1]}{$x}{'-'}){$p++;$x=-1;}else{last;}}
						if($x<0){
							while($i<$diff){$x = $a[1] + int(rand($a[2] - $a[1]));$st=int(rand()+0.5);if($st==1){print OUT "$a[0]\t",$x-($ext/2),"\t",$x-($ext/2)+$rl,"\tseq3\t0\t+\n";}else{print OUT "$a[0]\t",$x+($ext/2)-$rl,"\t",$x+($ext/2),"\tseq3\t0\t+\n";}$i++;}
						}
						$st=int(rand()+0.5);
						if(exists $bedF{$a[0]}{$a[1]}{$x}{$st}){if($st==1){$st=0;}else{$st=1;}}
						if($st==1){
						$bedF{$a[0]}{$a[1]}{$x}{$st}++;
						print OUT "$a[0]\t",$x-($ext/2),"\t",$x-($ext/2)+$rl,"\tseq3\t0\t+\n";
						$i++;
						}
						else{
						$bedF{$a[0]}{$a[1]}{$x}{$st}++;
						print OUT "$a[0]\t",$x-($ext/2),"\t",$x-($ext/2)+$rl,"\tseq3\t0\t+\n";
						$i++;
						}
					}
				}
			}
		}
		close NORM;
		close RAW;
		close OUT;
	}
	`rm $base.raw.bed`;
	print STDOUT "Done\n";
}
}

sub extractingGiveRegions{
print STDOUT "\n====Meta-profile plot module====\n";
$foldN=$reg;
if($reg==0 || $reg==""){print STDERR "\n\n***PLEASE PROVIDE REFERENCE FILE IN --plotRef OPTION***\n\n"; exit();}
$foldN=~s/.*\///;
#$foldN=$fileN;$foldN=~s/\..*//;
$fileN="Mid-$foldN";
$foldN=~s/\..*//;
$foldN="$opFolder/$foldN";
`mkdir -p $foldN`;
#`echo "Creating the reference using $reg file coordinates.....\$(date)"`;
foreach $hist (keys %fileList)
{
	foreach $samp (keys %{$fileList{$hist}})
	{
		$base="$opFolder/$hist/$samp/$hist-$samp";
		last;
	}
}


$time=localtime();
print STDOUT "\nChecking and preparing targets from $reg [$time]\n";
$opZ=`file $reg`;
if($opZ=~/compressed data/)
{
	$cmd="zcat $reg | $tooldir/bed_utils_light2 integrity $reg";
}
elsif($opZ=~/ASCII text/)
{
	$cmd="cat $reg | $tooldir/bed_utils_light2 integrity $reg";
}
else
{
	print STDERR "\n\n***PLEASE CHECK THE FILE. IT IS NEITHER COMPRESSED FILE NOR PLAIN TEXT FILE***\n\n"; exit();
}
		$intCode=`$cmd`;
		@intRes=split(",",$intCode);
		if($intRes[0]==1 || $intRes[0]>32)
		{
			print STDERR "\nError: Bed file is malformed in $reg file at this line $intRes[2]\n";
		}
		elsif($intRes[0]==2 || $intRes[0]==4)
		{
			print STDERR "\nError: Bed file is malformed in $reg file at this line $intRes[2].\n\tStrand should be in 6th column\n";
		}
		elsif($intRes[0]==8)
		{
			print STDERR "\nError: Bed file is malformed in $reg file at this line $intRes[2].\n\tSpace is used as delimiter, please use tab '\\t' as delimiter\n";
		}
		elsif($intRes[0]==16 && $intRes[0]==32)
		{
			print STDERR "\nError: Bed file is malformed in $reg file at this line $intRes[2].\n\tSpace is used as delimiter, please use tab '\\t' as delimiter\n\tStrand should be in 6th column\n";
		}
		if($intRes[0]){die();}
`awk '{mid=int((\$3+\$2)/2);print \$1"\t"mid"\t"mid"\t"\$4}' $reg | sort -k 1,1 -k 2,2n -k 3,3n -T $tmpdir | $tooldir/genetersect $reg $base.raw.bedgraph >$foldN/$fileN` if($plotT=~/target/i);
`awk '{if(\$6~/-/){print \$1"\t"\$3"\t"\$3"\t"\$4";"\$5":"\$2\$6\$3}else if(\$6~/\\+/){print \$1"\t"\$2"\t"\$2"\t"\$4";"\$5":"\$2\$6\$3}}' $reg | sort -k 1,1 -k 2,2n -k 3,3n -T $tmpdir | $tooldir/genetersect $reg $base.raw.bedgraph >$foldN/$fileN` if($plotT=~/tss/i || $plotT=~/genebody/i);

$time=localtime();
print STDOUT "\nReading target from $foldN/$fileN [$time]\n";
open(IN,"$foldN/$fileN") || die "Couldn't open the output reference file in specified folder $outF/$fileN\n";
$left=$left+$bin-($left) % $bin if($left % $bin);
$right=$right+$bin-($right) % $bin if($right % $bin);
while(<IN>)
{
	chomp();
	@a=split("\t",$_);
	@c=split(":",$a[3]);
	if($plotT=~/genebody/i)
	{
		@d=split("[+-]",$c[1]);
#		print "$d[0]\t$d[1]\n";
		if($c[1]=~/\+/){
			if($d[1]<$d[0]+($genePoints*$bin)){$st=$a[1]-1000;$end=$a[1]+(($genePoints*$bin)-1000);}
			else{$dif=($d[1]+500)-$d[0]; $st=$a[1]-1000;$end=$a[1]+$dif;$end=int($end+0.5);}
		}
		else{
			if($d[1]<$d[0]+($genePoints*$bin)){$st=$a[1]-(($genePoints*$bin)-1000);$end=$a[1]+1000;}
			else{$dif=$d[1]-($d[0]-500); $st=$a[1]-$dif;$st=int($st+0.5);$end=$a[1]+1000;}
		}
#	print "$st\t$end\t";
	$st=$st+$bin-($st) % $bin if(($st-1) % $bin);
	$end=$end+$bin-($end) % $bin if(($end-1) % $bin);
#	print "$st\t$end\t$a[3]\n";
	for($i=$st;$i<=$end;){$hash{$a[0]}{$i}{$a[3]}=1;$i=$i+$bin;}
	}
	else{
	if($c[1]=~/\+/){
		$st=$a[1]-$left+$bin;$end=$a[1]+$right;
	}
	else
	{
		$st=$a[1]-$right+$bin;$end=$a[1]+$left;
	}
	
	for($i=$st;$i<=$end;){$hash{$a[0]}{$i}{$a[3]}=1;$i=$i+$bin;}
	}
}
close IN;

foreach $k1 (sort keys %fileList)
{
	$time=localtime();
	print STDOUT "\n$k1: Building matrix [$time]\n";
	foreach $k2 (sort keys %{$fileList{$k1}})
	{
		$fileCur="$opFolder/$k1/$k2/$k1-$k2.$matTyp.bedgraph";
#		print "$fileCur.. "; system(date);
		$fileCurOP="$foldN/$k1-$k2.$matTyp.$plotT";
		$fileCurOPpng="$foldN/$k1-$k2.$matTyp.$plotT.png";
		open(CURF,$fileCur) || die "Couldn't open the file $fileCur\n";
		open(CURO,">$fileCurOP") || die "Couldn't open the file $fileCurOP\n";
		$mm=0;
		while(<CURF>)
		{
			chomp();
			if($_=~/^track type=bedGraph/){next;}
			@a=split("\t",$_);
#			@b=split("_",$a[3]);
			if(exists $hash{$a[0]}{$a[1]}){
				foreach $k (sort keys %{$hash{$a[0]}{$a[1]}}){
				$mmn=$a[3];
				$regions{$k}{$a[1]}=$mmn;
				}
			}
		}
		if($plotT=~/genebody/i){
			print CURO "RefID";
			for($i=1;$i<=$genePoints;$i++){print CURO "\tid$i";}
			print CURO "\n";
			foreach $mp (sort keys %regions)
			{
				print CURO "$mp";
#				print "$mp\t";
				undef @ck; undef @ckm; undef @ckm1; undef @val1; $loc=0;
				foreach $mpk (sort keys %{$regions{$mp}})
				{
					$mval=$regions{$mp}{$mpk};
					$idd="$k1-$k2";
					push(@val1,$mval);
				}
#				delete $regions{$mp};
				@c=split(":",$mp);
				if($c[1]!~/\+/){@val1=reverse(@val1);}
				$am=@val1;
				if($am>$genePoints){
					$dtp=int(1500/$bin);
					$cm=$am-$dtp;$cmd=int($cm/($genePoints-($dtp-1)));$loc=0;	#print "$mp\t$am\t$cm\t$cmd\n";
					for($mb=$dtp;($mb+$cmd-1)<=$am;)
					{
						$mc=$mb+$cmd-1;$loc++;if($loc==($genePoints-($dtp-1)-1)){$mc=$am-1;}
						@tmpd=@val1[$mb..$mc];$tmd=median(@tmpd); 
						undef @tmpd;
						push(@ckm1,$tmd);
						$mb=$mc+1;
					}
					@ckm=@val1[0..$dtp-1];push(@ckm,@ckm1);
				}
				elsif($am==$genePoints){
					@ckm=@val1;
				}
				else{
					print STDERR "Warning: Skipping $mp gene as it doesn't have enough length for $genePoints bins\n";
					next;
				}
				$ck1=join("\t",@ckm);$ck1=~s/\t$//g;
				print CURO "\t$ck1\n";
			}
			close CURF;
			close CURO;
		}
		else{
			print CURO "RefID";
			for($i=-($left);$i<$right;$i=$i+$bin){print CURO "\t$i";}
			print CURO "\n";
			foreach $mp (sort keys %regions)
			{
				print CURO "$mp";
				undef(@val1);
				foreach $mpk (sort keys %{$regions{$mp}})
				{
					$mval=$regions{$mp}{$mpk};
					$idd="$k1-$k2";
					push(@val1,$mval);
				}
				@c=split(":",$mp);
				if($c[1]!~/\+/){@val1=reverse(@val1);
				$mval=join("\t",@val1);
				print CURO "\t$mval\n";
				}
				else
				{
				$mval=join("\t",@val1);
				print CURO "\t$mval\n";
				}
			}
			close CURF;
			close CURO;
		}
	`Rscript $tooldir/MetaProfilePlot.R $fileCurOP $fileCurOPpng $plotT $matTyp Meta-profilePlot:$k1-$k2`;
	}
}
}


sub MAPlots {
print "\n====MA plot module====\n ";
foreach $k1 (sort keys %fileList)
{
	$time=localtime();
	print STDOUT "$k1: MA plotting [$time]\n";
	foreach $k2 (sort keys %{$fileList{$k1}})
	{
		@br=split("\t",$fileDen{$k1});
		if($denomList{$k1} eq $k2){next;}
		$fileOP="$opFolder/$k1"."_".$br[0]."_Vs_$k2.raw.png";
		$fileNum="$opFolder/$k1/$k2/$k1-$k2.raw.bedgraph";
		$denomFile="$opFolder/$k1/$br[0]/$k1-$br[0].raw.bedgraph";
		`Rscript $tooldir/MAPlot_v2.R $denomFile $fileNum $fileOP $k1 $br[0] $k2 BeforeNormalization`;
		$fileOPn="$opFolder/$k1"."_".$br[0]."_Vs_$k2.norm.png";
		$fileNumn="$opFolder/$k1/$k2/$k1-$k2.norm.bedgraph";
		$denomFilen="$opFolder/$k1/$br[0]/$k1-$br[0].norm.bedgraph";
		`Rscript $tooldir/MAPlot_v2.R $denomFilen $fileNumn $fileOPn $k1 $br[0] $k2 AfterNormalization`;
	}
}
}

sub zscore {
print STDOUT "\n====Z-score Normalization module====\n";
foreach $hist (sort keys %fileList)
{
	$time=localtime();
	print STDOUT "$hist: ZNorm running [$time]\n";
	foreach $samp (sort keys %{$fileList{$hist}})
	{
		$base="$opFolder/$hist/$samp/$hist-$samp";
		`Rscript $tooldir/zscore.R $base.norm.bedgraph $base.normZ.bedgraph $hist-$samp`;
	}
}
}

sub bedIntegrityCheck {
foreach $hist (sort keys %fileList)
{
	$time=localtime();
	print STDOUT "$hist: Checking bed file integrity [$time]\n";
	foreach $samp (sort keys %{$fileList{$hist}})
	{
		print STDOUT "\t$samp..";
		$fileCur=$fileList{$hist}{$samp};
		$opZ=`file $fileCur`;
		if($opZ=~/compressed data/)
		{
			$cmd="zcat $fileCur |";
	#		$txt="decompressing | counting..";
		}
		elsif($opZ=~/ASCII text/)
		{
			$cmd="cat $fileCur |";
	#		$txt="counting..";
		}
		else
		{
			print STDERR "\n\n***PLEASE CHECK THE FILE. IT IS NEITHER COMPRESSED FILE NOR PLAIN TEXT FILE***\n\n";
		}
		$cmd.="$tooldir/bed_utils_light2 integrity $fileCur";
		$intCode=`$cmd`;
		@intRes=split(",",$intCode);
#		print "$intCode\t$intRes[0]\n";
		if($intRes[0]==1 || $intRes[0]>32)
		{
			print STDERR "\nError: Bed file is malformed in $fileCur file at this line $intRes[2]\n";
		}
		elsif($intRes[0]==2 || $intRes[0]==4)
		{
			print STDERR "\nError: Bed file is malformed in $fileCur file at this line $intRes[2].\n\tStrand should be in 6th column\n";
		}
		elsif($intRes[0]==8)
		{
			print STDERR "\nError: Bed file is malformed in $fileCur file at this line $intRes[2].\n\tSpace is used as delimiter, please use tab '\\t' as delimiter\n";
		}
		elsif($intRes[0]==16 && $intRes[0]==32)
		{
			print STDERR "\nError: Bed file is malformed in $fileCur file at this line $intRes[2].\n\tSpace is used as delimiter, please use tab '\\t' as delimiter\n\tStrand should be in 6th column\n";
		}
		if($intRes[0]){die();}
		print STDOUT "Done\n";
	}
}
}

sub countReads {
print STDOUT "\n====Count Reads module====\n";
if($genome==0 || $genome==""){ print STDERR "\n\n***PLEASE PROVIDE REFERENCE CHROMSOME SIZE FILE IN --genomeInfo OPTION***\n\n"; exit();}
foreach $hist (sort keys %fileList)
{
	$time=localtime();
	print STDOUT "$hist: Counting reads [$time]\n";
	foreach $samp (sort keys %{$fileList{$hist}})
	{
	
		$fileCur=$fileList{$hist}{$samp};
		$id="$hist-$samp";$base="$opFolder/$hist/$samp/$id"; $histBase="$opFolder/$hist/$hist";
		print STDOUT "\t$samp..";
		$opZ=`file $fileCur`;
		if($opZ=~/compressed data/)
		{
			$opS=`zcat $fileCur | $tooldir/bed_utils_light2 issorted $fileCur`;
		}
		elsif($opZ=~/ASCII text/)
		{
			$opS=`cat $fileCur | $tooldir/bed_utils_light2 issorted $fileCur`;
		}
		else
		{
			print STDERR "\n\n***PLEASE CHECK THE FILE. IT IS NEITHER COMPRESSED FILE NOR PLAIN TEXT FILE***\n\n"; exit;
		}
		if($opZ=~/compressed data/)
		{
			if($opS==0)
			{
				$cmd="zcat $fileCur | sort -k 1,1 -k 2,2n -k 6,6 -T $tmpdir |";
				$txt="decompressing | sorting | counting..";
			}
			elsif($opS==1)
			{
				$cmd="zcat $fileCur |";
				$txt="decompressing | counting..";
			}
		}
		elsif($opZ=~/ASCII text/)
		{
			if($opS==0)
			{
				$cmd="sort -k 1,1 -k 2,2n -k 6,6 -T $tmpdir $fileCur |";
				$txt="sorting | counting..";
			}
			elsif($opS==1)
			{
				$cmd="cat $fileCur |";
				$txt="counting..";
			}
		}
		else
		{
			print STDERR "\n\n***PLEASE CHECK THE FILE. IT IS NEITHER COMPRESSED FILE NOR PLAIN TEXT FILE***\n\n";
		}
		print STDOUT " $txt";
		if($targetNorm){$reff="countbin";}else{$reff="wig";}
		$cmd="$cmd $tooldir/bed_utils_light2 $reff $fileCur $base.raw.bedgraph $genome $bin $extReads $hist-$samp";
		if($keepClonalReads==0)
		{
			if($noMiddlePos){$cmd="$cmd --noheader --noclonal";}
			else{$cmd="$cmd --noheader --noclonal --middlepos";}
		}
		else
		{
			if($noMiddlePos){$cmd="$cmd --noheader --clonal";}
			else{$cmd="$cmd --noheader --clonal --middlepos";}
		}
		$op=`$cmd`;
		print STDOUT "Done\n";
	}
}
}

sub quantileNormalization {
print STDOUT "\n====Quantile Normalization module====\n";
foreach $hist (sort keys %fileList)
{
	$histBase="$opFolder/$hist/$hist";
	$time=localtime();
	print STDOUT "$hist: QNorm running [$time]\n";
	undef @histFiles; undef @sampHeads; undef $sampNames; undef $histFileList;
	foreach $samp (sort keys %{$fileList{$hist}})
	{
		$base="$opFolder/$hist/$samp/$hist-$samp";
		push(@histFiles,"$base.raw.bedgraph");
		push(@sampHeads,"$hist-$samp");
	}
	$histFileList=join(" ",@histFiles);
	$sampNames=join("\t",@sampHeads);$sampNames="Chrom\tStart\tEnd\t$sampNames";
#	`awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") \$4 } END { for(i=1;i<=FNR;i++) print a[i] }' $histFileList >$histBase.raw.simple.bedgraph`;
	$cmd="cut -f 1-3";
	for($j=1;$j<=@histFiles;$j++){$pp=$j*4;$cmd.=",$pp";}
	$cmd="paste $histFileList | grep -v \"^track\" | $cmd | echo \"$sampNames\n\$(cat -)\" >$histBase.raw.bedgraph";
	`$cmd`;
	
	`Rscript $tooldir/quantileNorm_v2.R $histBase.raw.bedgraph $histBase.norm.bedgraph`;
	open(DT,"$histBase.norm.bedgraph") || die "Couldnt open the file $histBase.norm.bedgraph\n";
	while(<DT>)
	{
		chomp();
		$_=~s/\"//g;
		@aa=split("\t",$_);
		$lin="$aa[1]\t$aa[2]";
		if($_=~/^Chrom/)
		{
			$_=~s/\./-/g;
			for($i=3;$i<@aa;$i++)
			{
				@dt=split("\-",$aa[$i]);
				$aab="I$i";
				$fileSn="$opFolder/$hist/$dt[1]/$hist-$dt[1].norm.bedgraph";
				open($aab,">$fileSn") || die "Check1 Couldn't open the file $opFolder/$hg/$dtm/$hg-$dtm.norm.bedgraph\n";
				print $aab "track type=bedGraph name=$aa[$i]-QNorm\n";
			}
		next;
		}
		for($i=3;$i<@aa;$i++)
		{
			$ff=$ab[$i];
			$abb="I$i";
#			print $ff "$aa[0]\t$lin\t$aa[3]\t$aa[$i]\n";
			print $abb "$aa[0]\t$lin\t$aa[$i]\n";
#			$anom{$aa[0]}{$aa[3]}{$lin}{$ac[$i]}=$aa[$i];
#			print "$aa[0]\t$aa[3]\t$lin\t$ac[$i]\n" if($bnm<=3);
			$bnm++;
		}
	}
	$bnm=0;
	`rm $histBase.raw.bedgraph $histBase.norm.bedgraph` if($keepTemp==0);
}
}

sub usage {
$usageop="Tool:\tEpimetheus v1.0\nUsage:  perl $0 -c <config.file> -o <outputFolder> <[12345]|all> [options]\n\nGeneral options:\n--config | -c <file>\t\t-  Configuration File with filenames. Refer example. [Mandatory]\n--outputFolder | -o <string>\t-  Output folder name. [Mandatory]\n--storeTempFiles | -s \t\t-  Don't delete the temporary files.\n\nModules:\n--all\t-  Runs all the steps in one go\n--1\t-  Count Reads for all samples\n  --binSize | -b <integer>\t-  Window size; Should be multiples of hundred; [100]\n  --genomeInfo | -g <file>\t-  Chromosome sizes of corresponding genome build. [Mandatory for --all and --1]\n  --extReads | -e <int>\t\t-  Extend read length to mean fragment size [150]\n  --keepClonaReads | -k  \t-  Keep PCR duplicate reads in counting. By default it is omitted.\n  --noMiddlePos       \t-  By default middle position of the extended read will be used for read counts. Use this option to skip and use full reads for counting but normalized BED file cannot be generated with this option as one read will be counted in multiple bins.\n  --noIntegrityCheck       \t-  Skips the integrity check step (check the BED file integrity before counting step). Use this only when the files have been already subjected to integrity check by Epimetheus.\n--2\t-  Apply quantile normalization on results from step 1 and and z-score scaling on normalized values.\n  --normValueRound\t-  Round normalized value to nearest integer\n  --skipZscoreNorm\t-  Skip Z-score normalization; performs only quantile normalization\n--3\t-  Make MA plots between raw and normalized counts\n--4\t-  Make meta-profile plots on TSS centric or Genebody. Output will be generated in the the filename from --plotref inside given -o folder.\n  --plotType <genebody|tss|target>\t-  Plot whole genebody instead of flankings from TSS.\n  --plotDataFrom <raw|norm|normZ>\t-  Meta-profile plots should be plotted on raw or quantile normalized or Zscore normalized values.\n  --plotRef <file>\t\t\t-  A bed file. [Mandatory for --all and --4]\n  --left <integer>\t\t\t-  5' end flanking region size for plot\n  --right <integer>\t\t\t-  3' end flanking region size for plot\n--5\t-  Generate normalized BED file from normalized WIG file (Experimental), this step will give wrong output if the the previous runs were carried out with --noMiddlePos option\n  --readLenOut <integer>\t-  Length of the reads for normalized BED output [50]\n";
die "$usageop";
}

#Setting default values
$bin=100;
$keepTemp=0;
$help=0;
$matTyp="norm";
$left=1500;$right=1500;
$plotT="tss";
$extReads=150;
$opFolder=$config="";
$tmpdir="/tmp";
$readLenOut=50;

#Parsing options
GetOptions("1"=>\$countReads, "2"=>\$quantNorm, "3"=>\$makeMAPlots, "4"=>\$makeMetaprofilePlots, "5"=>\$wig2bed, "binSize=i"=>\$bin, "outputFolder=s"=>\$opFolder, "config=s"=>\$config, "genomeInfo=s"=>\$genome, "storeTempFiles"=>\$keepTemp, "extReads=i"=>\$extReads, "keepClonalReads"=>\$keepClonalReads, "help"=>\$help, "tmpdir=s"=>\$tmpdir, "plotType=s"=>\$plotT,,"plotRef=s"=>\$reg, "plotDataFrom=s"=>\$matTyp, "left=i"=>\$left, "right=i"=>\$right, "all"=>\$all, "noIntegrityCheck"=>\$integrityCheck, "pointsGenebody=i"=>\$genePoints, "noMiddlePos"=>\$noMiddlePos, "readLenOut=i"=>\$readLenOut, "skipZscorenorm"=>\$szsn);

if($all){
$countReads=1;$quantNorm=1;$wig2bed=1;$makeMAPlots=1;$makeMetaprofilePlots=1;
}

if($all==0 && $countReads==0 && $quantNorm==0 && $wig2bed==0 && $makeMAPlots==0 && $makeMetaprofilePlots==0)
{
	print STDERR "\n\n***PLEASE CHOOSE ANY MODULE TO RUN THE PIPELINE. FOR FIRST RUN, CHOOSE --ALL OPTION TO RUN THE WHOLE PIPELINE***\n\n";
	usage();
}

if($opFolder eq "" || $config eq "")
{
	print STDERR "\n\n***PLEASE PROVIDE MANDATORY OPTIONS --config [-c] AND --outputFolder [-o] CORRECTLY TO RUN THE PIPELINE***\n\n";
	usage();
}

if($makeMetaprofilePlots && $plotRef)
{
	print STDERR "\n\n***PLEASE PROVIDE --plotRef TO RUN THIS MODULE***\n\n";
	usage();
}

if($countReads && $genomeInfo)
{
	print STDERR "\n\n***PLEASE PROVIDE --genomeInfo <file> TO RUN THIS MODULE***\n\n";
	usage();
}

if($plotT!~/tss/i && $plotT!~/target/i && $plotT!~/genebody/i)
{
	print STDERR "\n\n***PLEASE PROVIDE tss OR genebody or target ONLY AS VALUES FOR --plotType***\n\n";
	usage();
	}
if($genePoints>0 && $plotT!~/genebody/i)
{
	print STDERR "\n\n***--pointsGenebody HAS TO BE USED ONLY WITH '--plotType genebody'***\n\n";
	usage();
}

if($readLenOut<35 || $readLenOut>200)
{
	print STDERR "\n\n***PLEASE PROVIDE AN OPTIMAL READ LENGTH [35-200]***\n\n";
	usage();
}

$genePoints=50;

if($help){usage();}
#Reading the configuration file to create hash with mark and sample
open(CFG, $config) || die "\nError: Couldn't open the config file $config; please check the usage or the file exists\n";
while(<CFG>)
{
	chomp();
	if($_=~/^#/){next;}
#	$_=~s/\s//g;
	($hist, $samp, $MAref, $fileName) = split ("\t",$_);
	$hist=~s/-/_/g;$samp=~s/-/_/g;
	if($MAref!="d" && $MAref!="n"){die "\nError: Please check in your configuration file $config; Only 'd' or 'n' is allowed in reference to be used for MA plot column\n";}
	if(exists $fileList{$hist}{$samp}){print STDERR "\nWarning: Please check in your configuration file $config; Sample $samp is mentioned twice for $hist; For each histone mark, sample names should be unique or only one will be processed\n";}
	$fileList{$hist}{$samp}=$fileName;
#Folder structure for each mark with each sample
	`mkdir -p $opFolder/$hist/$samp` if(!(-d "$opFolder/$hist/$samp"));
#To keep a list of reference files for MA plots to compare
	if($MAref=~/d/)
	{
		if(exists $denomList{$hist}){print STDERR "\nWarning: Please check in your configuration file $config; Mark $hist has two samples as reference for MA plots; For each histone mark, only one sample can be unique as reference for MA plots or choose one to one comparison option --MA1to1\n";}
		$fileDen{$hist}="$samp\t$fileName";
		$denomList{$hist}=$samp;
	}
}
close CFG;

if($countReads)
{
	bedIntegrityCheck() if($integrityCheck==0);
	countReads();
}
quantileNormalization() if($quantNorm);
zscore() if($quantNorm && $szsn==0);
MAPlots() if($makeMAPlots);
extractingGiveRegions() if($makeMetaprofilePlots);
wig2bed() if($wig2bed);
