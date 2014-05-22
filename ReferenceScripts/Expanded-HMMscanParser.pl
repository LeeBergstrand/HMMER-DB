#----------------------------------------------------------------------------------------------------------------------------------------
# Orginal Code.

# Yanbin Yin
# 08/18/2011
# hmmscan output parser
# Usage: sh hmmscan-parser.sh hmmscan-output-file

# 1. take hmmer3 output and generate the tabular output
# 2. sort on the 6th and 7th cols
# 3. remove overlapped/redundant hmm matches and keep the one with the lower e-values
# 4. calculate the covered fraction of hmm (make sure you have downloaded the "all.hmm.ps.len" file to the same directory of this perl script)
# 5. apply the E-value cutoff and the covered faction cutoff
cat $1 | perl -e 'while(<>){if(/^\/\//){$x=join("",@a);($q)=($x=~/^Query:\s+(\S+)/m);while($x=~/^>> (\S+.*?\n\n)/msg){$a=$&;@c=split(/\n/,$a);$c[0]=~s/>> //;for($i=3;$i<=$#c;$i++){@d=split(/\s+/,$c[$i]);print $q."\t".$c[0]."\t$d[6]\t$d[7]\t$d[8]\t$d[10]\t$d[11]\n" if $d[6]<1;}}@a=();}else{push(@a,$_);}}' \
	| sort -k 1,1 -k 6,6n -k 7,7n | uniq \| perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[0]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[2]<$c[2]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' \
        | uniq | perl -e 'open(IN,"all.hmm.ps.len");while(<IN>){chomp;@a=split;$b{$a[0]}=$a[1];}while(<>){chomp;@a=split;$r=($a[4]-$a[3])/$b{$a[1]};print $_."\t".$r."\n";}' \
	| perl -e 'while(<>){@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_ if $a[2]<1e-5;}else{print $_ if $a[2]<1e-3;}}' | awk '$NF>0.3'
''''
#--------------------------------------------------------------------------------------------------------------------------------------
# Expanded Perl

# cat $1 | perl -e 
while(<>)
{
	if(/^\/\//)
	{
		$x=join("",@a);
		($q)=($x=~/^Query:\s+(\S+)/m); # Grabs query name
		
		while($x=~/^>> (\S+.*?\n\n)/msg) # mutiline, "." matches newline, global
		{
			$a=$&;
			@c=split(/\n/,$a);
			$c[0]=~s/>> //;
			for($i=3;$i<=$#c;$i++)
			{
				@d=split(/\s+/,$c[$i]);
				print $q."\t".$c[0]."\t$d[6]\t$d[7]\t$d[8]\t$d[10]\t$d[11]\n" if $d[6]<1;
			}
		}
		@a=();
	}
	else
	{
		push(@a,$_);
	}
}

#| sort -k 1,1 -k 6,6n -k 7,7n | uniq \| perl -e
while(<>)
{
	chomp;
	@a=split;
	next if $a[-1]==$a[-2];
	push(@{$b{$a[0]}},$_); # Removes lines with same start and end
}

foreach(sort keys %b) 
{
        @a=@{$b{$_}};

        for($i = 0; $i < $#a; $i++) 
		{
                @b = split(/\t/, $a[$i]);   # Alignment 1 = top alignment
                @c = split(/\t/, $a[$i+1]); # Alignment 2 = bottom alignment
                $len1 = $b[-1] - $b[-2]; # Length one = aligment 1 end - aligment 1 start
                $len2 = $c[-1] - $c[-2]; # Length two = aligment 2 end - aligment 2 start
                $len3 = $b[-1] - $c[-2]; # Length one = aligment 1 end - aligment 2 start

                if($len3 > 0 and ($len3 / $len1 > 0.5 or $len3 / $len2 > 0.5)) # if alignments are overlaped and the overlap is greater than 50% the length of an alignment. 
                {
                        if($b[2] < $c[2]) # Checks E value. Removes the alignment with the lowest Evalue.
                        {
                                splice(@a, $i + 1, 1);
                        }
                        else
                        {
                                splice(@a, $i, 1);
                        }
                        $i = $i - 1;
                }
        }
        foreach(@a) 
		{
                print $_ . "\n";
        }
}

# | uniq | perl -e 
open(IN,"all.hmm.ps.len");
while(<IN>)
{
	chomp;
	@a=split;
	$b{$a[0]}=$a[1]; # creates hash of hmmName : hmmLength
}

while(<>)
{
	chomp;
	@a=split;
	$r=($a[4]-$a[3])/$b{$a[1]}; # $a[4] = hmm end $a[3] = hmm start ; $b{$a[1]} = result of the hash of the name of the hmm (hmm length).
	print $_."\t".$r."\n";
}
	
# | perl -e 
while(<>)
{
	@a=split(/\t/,$_);
	if(($a[-1]-$a[-2])>80)
	{
		print $_ if $a[2]<1e-5;
	}
	else
	{
		print $_ if $a[2]<1e-3;
	}
}
# awk '$NF>0.3' # Deletes alignment coverages less than 1/3

# ----------------------------------------------------------------------------------------------------------------------------------------
# Example Document to be parsed:

'''
# hmmscan :: search sequence(s) against a profile database
# HMMER 3.1b1 (May 2013); http://hmmer.org/
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             ../Blam.faa
# target HMM database:             Cyp125(GramPos).hmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (180 residues searched)
Target model(s):                           1  (408 nodes)
Passed MSV filter:                         0  (0); expected 0.0 (0.02)
Passed bias filter:                        0  (0); expected 0.0 (0.02)
Passed Vit filter:                         0  (0); expected 0.0 (0.001)
Passed Fwd filter:                         0  (0); expected 0.0 (1e-05)
Initial search space (Z):                  1  [actual number of targets]
Domain search space  (domZ):               0  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
Query:       YP_707205.1  [L=405]
Description: no_gene_name-hypothetical protein
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model    Description
    ------- ------ -----    ------- ------ -----   ---- --  -------- -----------

   [No hits detected that satisfy reporting thresholds]


Domain annotation for each model (and alignments):

   [No targets detected that satisfy reporting thresholds]

Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (405 residues searched)
Target model(s):                           1  (408 nodes)
Passed MSV filter:                         1  (1); expected 0.0 (0.02)
Passed bias filter:                        1  (1); expected 0.0 (0.02)
Passed Vit filter:                         0  (0); expected 0.0 (0.001)
Passed Fwd filter:                         0  (0); expected 0.0 (1e-05)
Initial search space (Z):                  1  [actual number of targets]
Domain search space  (domZ):               0  [number of targets reported over threshold]
# CPU time: 0.00u 0.00s 00:00:00.00 Elapsed: 00:00:00.00
# Mc/sec: inf
//
Query:       YP_703045.1  [L=110]
Description: no_gene_name-cytochrome P450
Scores for complete sequence (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Model           Description
    ------- ------ -----    ------- ------ -----   ---- --  --------        -----------
    1.3e-06   13.3   0.0    1.3e-06   13.3   0.0    1.0  1  Cyp125(GramPos)  


Domain annotation for each model (and alignments):
>> Cyp125(GramPos)  
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   13.3   0.0   1.3e-06   1.3e-06     142     188 ..      12      58 ..       1      77 [. 0.81

  Alignments for each domain:
  == domain 1  score: 13.3 bits;  conditional E-value: 1.3e-06
  Cyp125(GramPos) 142 aekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgedd 188
                      +++ ++d+v  +a+ +Pl+   + +G+p+  re+l+++ ++l+++  
      YP_703045.1  12 ETATTFDMVPALAAAFPLRVFPDAVGIPEVGRENLLSYGDHLFNAFG 58 
                      346789********************************999987655 PP
'''