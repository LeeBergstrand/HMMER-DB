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

#--------------------------------------------------------------------------------------------------------------------------------------
# Expanded Perl

# cat $1 | perl -e 
while(<>)
{
	if(/^\/\//)
	{
		$x=join("",@a);
		($q)=($x=~/^Query:\s+(\S+)/m); # Grabs query name
		
		while($x=~/^>> (\S+.*?\n\n)/msg) # 
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
	push(@{$b{$a[0]}},$_);
}
foreach(sort keys %b){@a=@{$b{$_}};

for($i=0;$i<$#a;$i++)
{
	@b=split(/\t/,$a[$i]);
	@c=split(/\t/,$a[$i+1]);
	$len1=$b[-1]-$b[-2];
	$len2=$c[-1]-$c[-2];
	$len3=$b[-1]-$c[-2];
	
	if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5))
	{
		if($b[2]<$c[2])
		{
			splice(@a,$i+1,1);
		}
		else
		{
			splice(@a,$i,1);
		}
		$i=$i-1;
	}
}
foreach(@a){print $_."\n";}} 

# | uniq | perl -e 
open(IN,"all.hmm.ps.len");
while(<IN>)
{
	chomp;
	@a=split;
	$b{$a[0]}=$a[1];}
	while(<>)
	{
		chomp;
		@a=split;
		$r=($a[4]-$a[3])/$b{$a[1]};
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

# ----------------------------------------------------------------------------------------------------------------------------------------
# Example Document to be parsed:

'''
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.1b1 (May 2013); http://hmmer.org/
# Copyright (C) 2013 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  ./TestHMMs/Cyp125(GramPos).hmm
# target sequence database:        -
# pfam-style tabular hit output:   /Users/lee/Dropbox/RandD/Repositories/HMM-Search-And-Extraction/TestGenomes/NC_008268rkdJTG.hmmer
# prefer accessions over names:    yes
# number of worker threads:        4
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Cyp125(GramPos)  [M=408]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence    Description
    ------- ------ -----    ------- ------ -----   ---- --  --------    -----------
   1.1e-231  767.3   0.0   1.3e-231  767.1   0.0    1.0  1  YP_704623.1  no_gene_name-cytochrome P450
   3.8e-133  442.8   0.0   4.2e-133  442.6   0.0    1.0  1  YP_704611.1  no_gene_name-cytochrome P450
   2.6e-130  433.4   0.1   3.2e-130  433.1   0.1    1.0  1  YP_702318.1  no_gene_name-cytochrome P450
   2.3e-112  374.3   0.0   2.7e-112  374.1   0.0    1.0  1  YP_702614.1  no_gene_name-cytochrome P450
    1.9e-95  318.6   0.0    2.4e-95  318.3   0.0    1.0  1  YP_703784.1  no_gene_name-cytochrome P450
    1.3e-75  253.3   0.0    1.1e-74  250.2   0.0    1.9  2  YP_704532.1  no_gene_name-cytochrome P450
    7.2e-66  221.2   0.0    9.2e-66  220.8   0.0    1.1  1  YP_704613.1  no_gene_name-cytochrome P450
      9e-60  201.1   0.0    1.2e-59  200.7   0.0    1.0  1  YP_705656.1  no_gene_name-cytochrome P450
    1.1e-57  194.2   0.0    1.9e-57  193.4   0.0    1.3  1  YP_700387.1  no_gene_name-cytochrome P450
    6.4e-56  188.4   0.0    9.3e-56  187.9   0.0    1.1  1  YP_702473.1  no_gene_name-cytochrome P450
    7.8e-54  181.5   0.0    1.1e-53  181.0   0.0    1.3  1  YP_700417.1  no_gene_name-cytochrome P450
    8.5e-54  181.4   0.0    1.1e-53  181.1   0.0    1.2  1  YP_702911.1  no_gene_name-cytochrome P450
    6.7e-53  178.5   0.0    8.1e-53  178.2   0.0    1.0  1  YP_702567.1  no_gene_name-cytochrome P450
    9.4e-52  174.7   0.0    1.6e-51  173.9   0.0    1.3  1  YP_700371.1  no_gene_name-cytochrome P450
    6.4e-30  102.8   0.0    9.1e-30  102.3   0.0    1.1  1  YP_702345.1  no_gene_name-cytochrome P450
    1.2e-26   92.0   0.0    1.9e-26   91.4   0.0    1.3  1  YP_703834.1  no_gene_name-cytochrome P450
    1.5e-11   42.3   0.0    1.7e-11   42.1   0.0    1.0  1  YP_702742.1  no_gene_name-cytochrome P450
    6.9e-10   36.8   0.2    5.9e-08   30.4   0.0    2.8  2  YP_705149.1  no_gene_name-cytochrome P450
    1.2e-09   36.0   0.1    0.00034   18.1   0.0    3.0  3  YP_703037.1  no_gene_name-cytochrome P450
    2.5e-09   35.0   0.0    7.9e-07   26.7   0.0    2.1  2  YP_704615.1  no_gene_name-cytochrome P450
    0.00012   19.6   0.0     0.0041   14.5   0.0    2.4  2  YP_704571.1  no_gene_name-cytochrome P450
     0.0093   13.3   0.0     0.0096   13.3   0.0    1.0  1  YP_703045.1  no_gene_name-cytochrome P450
  ------ inclusion threshold ------
      0.058   10.7   0.0      0.064   10.6   0.0    1.0  1  YP_702741.1  no_gene_name-hypothetical protein


Domain annotation for each sequence (and alignments):
>> YP_704623.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  767.1   0.0  4.1e-234  1.3e-231       2     408 .]      55     460 ..      54     460 .. 1.00

  Alignments for each domain:
  == domain 1  score: 767.1 bits;  conditional E-value: 4.1e-234
  Cyp125(GramPos)   2 aapalpegfdftdpdllaerlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielq 95 
                      a+p+lpegfdftdpd++aer+P++e+aelrktaP+wWn+qp e+ggf+d+GyWvv+k +dvkevsrrsdvfs++entaivrf+ddi re+ie+q
      YP_704623.1  55 AQPNLPEGFDFTDPDVYAERIPYQEFAELRKTAPIWWNPQPPEIGGFHDDGYWVVSKLEDVKEVSRRSDVFSTHENTAIVRFADDIPRENIEMQ 148
                      689******************************************************************************************* PP

  Cyp125(GramPos)  96 klvllnkdapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddp 189
                      +++l+nkdapeht+lrk+vsr+ftPra+++lreel+era+kivk+aae+g+gdfv+qva+elPlqaiaellGvpqedr k+fdwsn+++g+ddp
      YP_704623.1 149 RFILINKDAPEHTKLRKLVSRGFTPRAINSLREELTERAEKIVKEAAESGAGDFVTQVACELPLQAIAELLGVPQEDRLKVFDWSNQMTGYDDP 242
                      ********************************************************************************************** PP

  Cyp125(GramPos) 190 efaeedaaaasaellayamklaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrerpe 283
                      e+ + d++aas+e+l+ya+++a+erkk+Paddivt+l+eadidg++ls++efgffvillavaGnettrnaithGm+afld+pdqWelyk+erp+
      YP_704623.1 243 EL-DIDPQAASMEILGYAYQMADERKKCPADDIVTTLIEADIDGNELSPEEFGFFVILLAVAGNETTRNAITHGMMAFLDHPDQWELYKKERPK 335
                      *9.******************************************************************************************* PP

  Cyp125(GramPos) 284 taadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifnai 377
                      t+adeivrwatPv+sfqrtaledtelggv+ikkgqrvv++y sanfde+ f++Pe+fdi+r++nph+gfGgtGah+c+Ganlarleidlifnai
      YP_704623.1 336 TTADEIVRWATPVNSFQRTALEDTELGGVQIKKGQRVVMLYGSANFDEDAFENPEKFDIMRENNPHVGFGGTGAHFCLGANLARLEIDLIFNAI 429
                      ********************************************************************************************** PP

  Cyp125(GramPos) 378 adklPdlkklgeperlrsgwlngikelqvdy 408
                      ad+lPd++klg+p+rlrsgwlngike+qvdy
      YP_704623.1 430 ADHLPDISKLGDPRRLRSGWLNGIKEFQVDY 460
                      *****************************98 PP

>> YP_704611.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  442.6   0.0  1.4e-135  4.2e-133       8     407 ..       4     398 ..       1     399 [. 0.95

  Alignments for each domain:
  == domain 1  score: 442.6 bits;  conditional E-value: 1.4e-135
  Cyp125(GramPos)   8 egfdftdpdllaerlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvlln 101
                      + +d+  pd+++e++P++ + +lr ++Pv+W+++++++      G+W+vt+h+dv  vsr+s +fss+  t+ +   dd ++++   q ++l+n
      YP_704611.1   4 SRIDLKCPDVYTEGVPYAFFDHLRSSEPVYWQPEENGT------GFWAVTRHADVVAVSRDSVTFSSAVGTTQI---DDFDEQTRAKQAAMLVN 88 
                      5689999************************9998875......9**********************9877643...44444555567789*** PP

  Cyp125(GramPos) 102 kdapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefa..e 193
                       d+p+htrlr++vsr+ftPr+v+ l+ ++++ +++iv++  e+   dfv + a+ lPl+ ia llG p  d ++l+dwsn+++g ddpe+   +
      YP_704611.1  89 LDPPDHTRLRQLVSRGFTPRTVKVLEGHIRDICTRIVDRLLEARDVDFVPEAAAPLPLEVIAALLGAPPGDVDRLYDWSNRMIGFDDPEYGttQ 182
                      *****************************************************************************************97557 PP

  Cyp125(GramPos) 194 edaaaasaellayamklaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrer..peta 285
                      +d++ a+ae++ ya +la++r+ +P ddivtklv+ d +g++l++ ef +f +ll+vaGnettrna++ Gm+af+d+pdqW+  + e   +++a
      YP_704611.1 183 ADGELAAAEIFLYANELAAQRRIEPRDDIVTKLVQPDENGDTLTEMEFNMFFVLLVVAGNETTRNATAGGMQAFIDHPDQWRRLQSEPdlASSA 276
                      789*******************************************************************************999875226899 PP

  Cyp125(GramPos) 286 adeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifnaiad 379
                      ++ei+rw tPv+ f+rta++d  +g++ ++ g++vv++y san de vfd+P++fdi r+pnp ++fGg+G h+c+Ga+larle++++f+ +a 
      YP_704611.1 277 VEEILRWVTPVMDFRRTATRDAYIGDQLVRAGDKVVIYYPSANRDEAVFDEPYRFDIGRSPNPQIAFGGSGVHFCLGAHLARLELRILFETLAA 370
                      ********************************************************************************************** PP

  Cyp125(GramPos) 380 klPdlkklgeperlrsgwlngikelqvd 407
                      ++  ++++g   rlrs ++ gik+++v+
      YP_704611.1 371 RIDRVESTGPVARLRSNFISGIKTMPVR 398
                      ************************9986 PP

>> YP_702318.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  433.1   0.1    1e-132  3.2e-130       9     407 ..       5     401 ..       2     402 .. 0.94

  Alignments for each domain:
  == domain 1  score: 433.1 bits;  conditional E-value: 1e-132
  Cyp125(GramPos)   9 gfdftdpdllaerlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvllnk 102
                      ++d+ +pd +a+++P+e +a lr+ aPv+ +  +++      + +W vt+h+d+ +v r+++++ss++ +    + dd++ +++  q+l++ln+
      YP_702318.1   5 DVDLYNPDTFAKGVPHEMFAVLRREAPVYRHLDERG------DPFWCVTRHADIVTVNRDAETYSSWRGA---TYIDDLSPDDLAGQQLMMLNM 89 
                      789*************************97766655......68***********************986...5889***************** PP

  Cyp125(GramPos) 103 dapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefaeed. 195
                      d+p+ht lrkivs++ftPr + +l+e l++ra++iv+a  e+g+ dfv +va+elPlqaia++lGvpqedr+ +fd +n+++g++dpef  ed 
      YP_702318.1  90 DPPDHTALRKIVSKGFTPRRIGQLHEILARRATTIVDAVIERGECDFVVDVASELPLQAIADFLGVPQEDRKLIFDLTNQMIGSSDPEFHLEDg 183
                      *****************************************************************************************96664 PP

  Cyp125(GramPos) 196 .aaaasaellayamklaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykr..erpetaa 286
                       ++aa+a+++ay++++ e+rkk+P ddi t l++ad++gekl + +f +f +llavaGnettrnai+h  la++++p+  +   +  ++ ++ +
      YP_702318.1 184 qERAAAAQMFAYSLEMFEDRKKHPRDDIATALIQADVHGEKLGELDFNMFFMLLAVAGNETTRNAISHTQLALMEHPEERRKVLEdpSKLDALI 277
                      4678999**********************************************************************98876554115667889 PP

  Cyp125(GramPos) 287 deivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdp...nphlgfGgtGahyciGanlarleidlifnai 377
                      +e +rwatPv+ f+rta++dt l +v+i +g+rvv+ + s n de vfddP+tfdi r       h++fGg G h+c+Ganlar e+++++ +i
      YP_702318.1 278 EEGLRWATPVMQFRRTATTDTVLHDVEISEGDRVVIWHMSGNRDEAVFDDPYTFDIDRPTghySQHIAFGGGGHHFCLGANLARAEMKVMLSEI 371
                      *********************************************************96522268***************************** PP

  Cyp125(GramPos) 378 adklPdlkklgeperlrsgwlngikelqvd 407
                        ++Pd+++++ ++rlrs ++ng+k+++v+
      YP_702318.1 372 LRRMPDMEQTAPAQRLRSNFINGLKHMPVT 401
                      ***************************997 PP

>> YP_702614.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  374.1   0.0  8.6e-115  2.7e-112      10     407 ..       8     404 ..       2     405 .. 0.95

  Alignments for each domain:
  == domain 1  score: 374.1 bits;  conditional E-value: 8.6e-115
  Cyp125(GramPos)  10 fdftdpdllaerlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvllnkd 103
                       d+  pd + +++P++ +aelr+ +Pv W+e+p++ g     G+W+vt+++dv +vs+++dvfss++ ++ +r   d + +++ + ++++ln d
      YP_702614.1   8 TDLGSPDAYIHGVPHQVFAELRRHEPVAWIEEPAGEGFAGGPGFWAVTRYDDVMTVSKKPDVFSSHKGASFLR---DQSPQDLAALQQMMLNLD 98 
                      58899******************************9877889************************9988766...5567788888999***** PP

  Cyp125(GramPos) 104 apehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefaeed.a 196
                      +p+h+++r ivs++ftP+ v+ + ++++++ar+iv+a  ++g+ d+ve+v++e+Pl+ +a++lGvp+edr+ l+dw+n++vg ddp++  ++  
      YP_702614.1  99 PPDHSQMRSIVSKVFTPKMVRGMFDSIADHARAIVDALPDDGEIDLVEHVSAEMPLRVLADVLGVPSEDRHLLYDWTNRMVGLDDPSYGGREaF 192
                      ****************************************************************************************988824 PP

  Cyp125(GramPos) 197 aaasaellayamklaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrerp..etaade 288
                       +a +e++ y+ + + +++++P+ d+ + +v+a++dg++ls++e+  f  ll++aGnettrn +t  +l++ ++p   e    + +    a++e
      YP_702614.1 193 LSAFIEMFEYSAAQTRAKRTEPGSDVWSLIVNAEVDGTQLSPEELNRFFQLLVIAGNETTRNLLTGAILTLGEHPGEREKLADDPAllPNAIEE 286
                      56779************************************************************************99998765422589*** PP

  Cyp125(GramPos) 289 ivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifnaiadklP 382
                      ++r+ +Pv+ f+rt ++ telgg+++++gq+vv+fy san de  fddP+tf i r    hl+fG  G+h+c+G++larle +++++a+  ++P
      YP_702614.1 287 VLRFHSPVIQFRRTVTRSTELGGKQMHEGQKVVIFYVSANRDEAQFDDPDTFRIDRGAANHLAFG-AGTHFCLGNSLARLEAKVLLEALFTRFP 379
                      ****************************************************************8.69************************** PP

  Cyp125(GramPos) 383 dlkklgeperlrsgwlngikelqvd 407
                      + + +g p+r+rs ++ gik+l+v+
      YP_702614.1 380 HWQVTGPPDRFRSNFIHGIKKLPVH 404
                      **********************996 PP

>> YP_703784.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  318.3   0.0   7.6e-98   2.4e-95      26     405 ..      28     402 ..      19     405 .. 0.93

  Alignments for each domain:
  == domain 1  score: 318.3 bits;  conditional E-value: 7.6e-98
  Cyp125(GramPos)  26 elaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvllnkdapehtrlrkivsrlft 119
                       + elr+  P+ +++++++       G+W v  ++dv  vsr+++vfss++  +i    dd+ +e +e + + ++ +d+p+h rlr iv+ +ft
      YP_703784.1  28 MFEELRQRSPIEFHAETAG------PGFWSVVGYDDVVAVSRNPEVFSSAKGFTI----DDVPAEILEFA-MSMIAMDDPRHRRLRGIVQSAFT 110
                      5789999999999877655......59**********************998887....67777778876.5799******************* PP

  Cyp125(GramPos) 120 PravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelv.geddpefaeed.aaaasaellayamkla 211
                        +v+ + e+++  a +i+++  ++ ++dfv +va++lP+q i+el+G+p++dr  ++  + ++v g  d+ef ++  +aa   ++ +ya++l 
      YP_703784.1 111 AASVRGIAERVRGTAAQIAADLPRDAEFDFVPDVATRLPVQVICELIGIPENDRPMILAAAADVVaGGGDQEFVTAPgGAAGLGQIYGYALELG 204
                      *********************************************************9988777626789999877636677789********* PP

  Cyp125(GramPos) 212 eerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykre...rpetaadeivrwatPvtsfqrt 302
                      e+rk+ P++d+ ++lv++ +dge+l+++efg fvill  aG ettr a++  +l + ++pdq el   +     + ++de++r+a+Pv  ++rt
      YP_703784.1 205 ERRKADPGEDLTSRLVSTMVDGEALAPTEFGSFVILLITAGFETTRQALAWALLLLSEHPDQKELLLADfagHIDNTIDEVIRFASPVPYMRRT 298
                      *****************************************************************987622256789***************** PP

  Cyp125(GramPos) 303 aledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilr.dpnphlgfGgtGahyciGanlarleidlifnaiadklPdlkklgeperlrs 395
                      a+ dtelgg +i  g++vv+ y san+de vfd P +fdi+r +   hlgfG++  h+c+G nlarle++++++++  + P l+ +gepe l s
      YP_703784.1 299 ATVDTELGGAHIAAGDKVVMWYLSANHDERVFDGPGRFDITRaNAGKHLGFGAKDIHHCLGVNLARLELRVMLEELFTAHPTLHAVGEPELLLS 392
                      ******************************************44568*********************************************** PP

  Cyp125(GramPos) 396 gwlngikelq 405
                      +++ g+ +l+
      YP_703784.1 393 AFISGVAALP 402
                      ******9876 PP

>> YP_704532.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?    1.0   0.0      0.16        50       7      43 ..      10      48 ..       4      55 .. 0.80
   2 !  250.2   0.0   3.6e-77   1.1e-74      99     407 ..      82     396 ..      79     397 .] 0.92

  Alignments for each domain:
  == domain 1  score: 1.0 bits;  conditional E-value: 0.16
  Cyp125(GramPos)  7 pegfdftdpdllaerl..PveelaelrktaPvwWneqpe 43
                      ++fd++d   ++  l  P + +a +r  +Pv+  +   
      YP_704532.1 10 RPDFDLIDGRFYSGELgdPRQAYAWMRAHEPVYRADGIL 48
                     579*********98877799************9776655 PP

  == domain 2  score: 250.2 bits;  conditional E-value: 3.6e-77
  Cyp125(GramPos)  99 llnkdapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefa 192
                      ++++d+p+h + r++v+ +ft + vea   +++e + ++++aa +k++ dfv+++a+ lP+  + ++lG+  e+r+  ++ws++l+++  ++ +
      YP_704532.1  82 MIDMDDPQHLQHRRLVNAGFTRKQVEAKIGRIREICDHLIDAACDKDEVDFVRDLAAPLPMAVVGDMLGMRPEERDTFLQWSDDLMNALGSNAT 175
                      89************************************************************************************99988887 PP

  Cyp125(GramPos) 193 eed.aaaasaell..ayamklaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykre..r 281
                       e+ +a+a+a l    ++++  eer++nP+dd+ + lv+++idg++lsd ++    +l+ + G+ettr  ++ Gm  ++ +pdq +   ++   
      YP_704532.1 176 PEElQAQANAYLAfnEFTLRTIEERRQNPTDDLTSILVHSEIDGHRLSDPDIVGETLLILIGGDETTRHVLSGGMEQLMRHPDQHQRLVNDpdG 269
                      5554555555544226788899******************************9999999**************************988765114 PP

  Cyp125(GramPos) 282 petaadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifn 375
                        +a++e++rw++P++ + rt ++d+e+ g+ +++g++++l+++sanfd+ vfd+Pe fdi r+pnph++fG  G+h+c+G+++arle +++f+
      YP_704532.1 270 IPAAVEEMLRWSSPIKNMCRTVTRDVEFFGTDLRQGEKMMLLFESANFDDAVFDSPEVFDIERSPNPHVAFG-FGTHFCLGNQFARLEAKIMFE 362
                      5689*******************************************************************7.79******************* PP

  Cyp125(GramPos) 376 aiadklPdlkklge.p.erlrsgwlngikelqvd 407
                      ++  +lPd+  +++ p +r  + ++ g++e++v+
      YP_704532.1 363 QLLSRLPDMVLVDDgPlRRRPANFVSGLEEMPVK 396
                      **********9987442455578999*9999986 PP

>> YP_704613.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  220.8   0.0   2.9e-68   9.2e-66      27     404 ..      31     399 ..      14     403 .. 0.89

  Alignments for each domain:
  == domain 1  score: 220.8 bits;  conditional E-value: 2.9e-68
  Cyp125(GramPos)  27 laelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvllnkdapehtrlrkivsrlftP 120
                      + +lr+ aP+++ne+           +W++++h+dv    r+++  ss++  ++       +      + + +l  d+p+h r+r++vsr+ftP
      YP_704613.1  31 YRRLREEAPLYYNEE---------LDFWALSRHEDVVAAFRDNQRLSSANGVSLDPA----AYGPHAHKTMSFLALDDPRHMRMRQLVSRGFTP 111
                      555555555555543.........359*******************99887766443....344445566789999****************** PP

  Cyp125(GramPos) 121 ravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefa.eedaaaasaellayamklaee 213
                      r v++l+ ++ + +r+  + a ++g++d++ + a +lP+  i+el+Gvp++dr +l + ++ +v +++  +   + a++as  l++y  ++ +e
      YP_704613.1 112 RRVNELEGRILDLTRQYLDPALAAGEFDWIGEFAGKLPMDVISELMGVPEADRVELRRLADLVVHREEGVLDvPHAAMEASLYLVGYYADMLAE 205
                      *****************************************************************99998864566899*************** PP

  Cyp125(GramPos) 214 rkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqW.elyk.rerpetaadeivrwatPvtsfqrtale 305
                      r+++P++d+ + l+ea+idg++l+dde+  f+ l++vaGnett   + + +   + +p q  +++   +r    ++e +r+ t  + + rta  
      YP_704613.1 206 RRRKPTEDLTSALLEAEIDGDRLTDDEIIAFMFLMVVAGNETTTKLLGNALYWGFRYPGQArNVFEdASRVPEWVEETLRFDTSSQMVARTASV 299
                      ***************************************************99999999862578625788999******************** PP

  Cyp125(GramPos) 306 dtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifnaiadklPdlkkl.geperlrsgwl 398
                      d ++ +++i  g++v ++  san d  vf+d +++ i rd    l+  g G+h+c+Ga++arle ++ + ++  ++ d     +++ r++s+ +
      YP_704613.1 300 DLDFHDRTIPAGDKVLILIGSANRDSAVFEDADEYRIGRDTSNKLASFGGGTHFCLGAHMARLEARIALTELVSRIDDYDIDeANSVRVHSTNV 393
                      *****************************************999988888**************************997654145678888888 PP

  Cyp125(GramPos) 399 ngikel 404
                       g  +l
      YP_704613.1 394 RGFATL 399
                      887666 PP

>> YP_705656.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  200.7   0.0   3.8e-62   1.2e-59      21     385 ..      28     388 ..      15     396 .. 0.88

  Alignments for each domain:
  == domain 1  score: 200.7 bits;  conditional E-value: 3.8e-62
  Cyp125(GramPos)  21 rlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqk.lvllnkdapehtrlrki 113
                      r P + +a lr  +Pv+ +   e++ g   + yWv+++h+dv   +r++++fss++    + + +    e+i lq    l++ d+p+ht +r++
      YP_705656.1  28 RAPWQMYAGLRDHDPVHHVVP-EDAPG---NDYWVLSRHADVYAAARDPETFSSAAG-LTTTYGE---LEKIGLQDnPPLVMLDPPDHTAFRRL 113
                      7899************98765.44433...56*********************9864.4444544...455555541356777*********** PP

  Cyp125(GramPos) 114 vsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefaeedaaaasaellaya 207
                      vs++ftPr v+a++ ++++   + +++ +e+g+gd v+++   lP   +a+ lGvp+edr++   w++++v+++ +   + +a+ a +el++y 
      YP_705656.1 114 VSKGFTPRQVTAVEPNVRAFVVERIERLREAGAGDVVKELFKPLPSMVVAHYLGVPDEDRAQFDGWTEAIVAANAQGD-TLQATGAVTELMGYF 206
                      ***********************9***********************************************9987775.89999********** PP

  Cyp125(GramPos) 208 mklaeerkknPaddivtklveadi..dgeklsddefgffvillavaGnettrnaithGmlafldnpdqWe.lykr.erpetaadeivrwatPvt 297
                        l e+r++ P+dd +++lv++ +  dg+      +  f   ++  Gn+tt   +   ++ ++++ dq + l  + e  + a++e++r ++Pv+
      YP_705656.1 207 TGLIERRRTDPGDDTISHLVAGGMgaDGDITGLLSILGFAFTMVTGGNDTTTGMLGGAVQLLTEHRDQRRdLIEHpELVRDAVEELLRLTAPVQ 300
                      ********************98654377765555655566677889**********************8625443267789************* PP

  Cyp125(GramPos) 298 sfqrtaledtelggvkikkgqrvvlfyrsanfdeevfd.dPetfdilrdpnphlgfGgtGahyciGanlarleidlifnaiadklPdlk 385
                       + rt ++d+e+ g+ + +g +v l+y san d+  f  d e +d+ r+p   l+f ++Gah+c+Ga  ar++ ++ ++++  + Pd++
      YP_705656.1 301 GLARTVTRDVEIEGTVVPEGRKVLLLYGSANRDPRRFGpDAEVLDVRRSPKQILTF-SHGAHHCLGAAAARMQARIALEELLARCPDFT 388
                      *************************************7478899************.58***************************986 PP

>> YP_700387.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  193.4   0.0   6.1e-60   1.9e-57       5     407 ..      11     406 ..       7     407 .. 0.87

  Alignments for each domain:
  == domain 1  score: 193.4 bits;  conditional E-value: 6.1e-60
  Cyp125(GramPos)   5 alpegfdftdpdlla.erlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqkl 97 
                      +l ++fd+ dp+l     +  e +a lrk+ Pv + ++         +G+Wvvt++k+  +v +++++fss  n  +             + k+
      YP_700387.1  11 HLITDFDVYDPSLAVpADVFQERVAALRKQGPVLYSPH--------HGGHWVVTRYKEALQVLQDPETFSSFPNNLLNA----------AQGKF 86 
                      56789*******98636789999***********9875........5799*********************99987644..........35789 PP

  Cyp125(GramPos)  98 vllnkdapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelv....ged 187
                      + l+ d+peh+  r+ ++ lf P+ ++al+ e+++  ++++++ a++g+ +f+++ a elP +    l+G p +d  +  +w++  +    g++
      YP_700387.1  87 LPLELDPPEHSYYRQALQPLFSPKQMKALEPEIRKIITELIDQFADRGECEFISEFAHELPTRIFLALMGWPLSDAPQFTEWTDITLqgipGAS 180
                      99*********************************************************************************85441111344 PP

  Cyp125(GramPos) 188 dpefaeedaaaasaellayamkl.aeerk.knPaddivtklvea..didge..klsddefgffvillavaGnettrnaithGmlafldnpdqWe 275
                      ++e  +e  a+a+ e+  y  k+ a  r+ +  ++ + ++++++  didg+  +l+d+e++ +  ll vaG  t + a++ G++ + +n +q e
      YP_700387.1 181 EAES-AEARAKAAGEIYEYFGKVvARVRSgEDSSESLTAQIINTplDIDGTprSLTDEELSRMFFLLLVAGLHTVQGALAWGLIHLSHNLEQRE 273
                      4444.4446778899999976651666652455677888888875577885448**************************************** PP

  Cyp125(GramPos) 276 lykrerp..etaadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlar 367
                          +     +a++ei+r  t  +   r a++d e+ggv+ik g+++ ++ +san d++ fd+Pe+++i r+pn h+gfG  G+h c+G++lar
      YP_700387.1 274 AIIDDPTliPSAVEEILRIETATS-SGRRATRDAEIGGVSIKAGDQLLVVMTSANRDNDEFDSPEELQIERHPNRHIGFG-AGPHRCLGSHLAR 365
                      987765422689*******99865.57889*************************************************8.69*********** PP

  Cyp125(GramPos) 368 leidlifnaiadklPdlkkl.geperlrsgwlngikelqvd 407
                      le++l +++i  ++Pd + + ++p  ++s+ + g ++l+++
      YP_700387.1 366 LELRLAMEEIHRRIPDYALVpDNPAVFHSSQVRGCEKLPIT 406
                      *****************9984679999***99999998876 PP

>> YP_702473.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  187.9   0.0     3e-58   9.3e-56      47     394 ..      32     387 ..      15     401 .. 0.87

  Alignments for each domain:
  == domain 1  score: 187.9 bits;  conditional E-value: 3e-58
  Cyp125(GramPos)  47 gfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvllnkdapehtrlrkivsrlftP....ravealreeleerark 136
                      + +d+G ++v ++ dvk + +++ + s  +n      +    +++       +l+ d+peh rlr++ +r+f P      ++++r el++  +k
      YP_702473.1  32 APQDNGSFLVGRYYDVKALLHDPRISSDLHNRG----PGLPGADRAGDGPESFLKLDPPEHDRLRRLTTRQFGPphspDRIDSMRGELAALVTK 121
                      567889999999999999888888877766643....33334555556667899******************75222268************** PP

  Cyp125(GramPos) 137 ivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefaeed......aaaasaellayamklaeerkknPaddivt 224
                      +v+    k++ d v++ a+ +P+  i+ llGvp+e+  + + w+++lv++ dp  a +d      a++a++e+ ++   l e+r++ P+dd+++
      YP_702473.1 122 LVDGLDGKERIDIVDDFAYPFPVTVICRLLGVPHEEEPRFHAWADALVAAIDPRRAGADneqvkvAQKAQMEMAMFMAGLIEDRRRDPGDDMLS 215
                      *****************************************************9985552233334568889999999999************* PP

  Cyp125(GramPos) 225 klveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrer..petaadeivrwatPvtsf.qrtaledtelggvkik 315
                       l +    + +l+  e+    ill +aG ett n it+Gml++l npd  e  + e      a++e++r+  Pv+ + qrt + d+e+ggv++ 
      YP_702473.1 216 GLANDSGPDGQLNTVELVTNSILLFIAGHETTVNLITNGMLTLLRNPDVLERLRAEPelMPQAVEELLRFEPPVHLLgQRTPIVDIEIGGVTVP 309
                      ***98887778*****************************************99875114579***********9877**************** PP

  Cyp125(GramPos) 316 kgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifnaiadklPdlkklgeperlr 394
                      kg  +vl  +san d+e f d ++f   r +n h+g+G +G h c Ga larle ++ + ++  +lP+ + +++p   r
      YP_702473.1 310 KGSPLVLALASANRDPERFVDADKFVPDRADNQHVGLG-SGIHSCFGAPLARLEGQMALAELIRRLPNPRLVEDPPPYR 387
                      ************************************97.6****************************99999986555 PP

>> YP_700417.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  181.0   0.0   3.6e-56   1.1e-53       8     404 ..      12     401 ..       1     405 [. 0.83

  Alignments for each domain:
  == domain 1  score: 181.0 bits;  conditional E-value: 3.6e-56
  Cyp125(GramPos)   8 egfdftdpdlla.erlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqklvll 100
                       +fd+ d  l + e +  e  a lr   Pv++ +          +G+W+vt+++++++v r++++fss+ n  +             + k++ +
      YP_700417.1  12 VDFDVYDQTLAMpEDVFQERAAALRAIGPVVYSKA--------HGGHWIVTRYEEIHQVLRDPETFSSYPNNLVN----------AGQGKFIPI 87 
                      4899999888653678888899999999*999875........5799*********************9998774..........356899*** PP

  Cyp125(GramPos) 101 nkdapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefa.. 192
                      + d+peht  r+ ++ lf P+ +++l+ ++++  ++++++ a++g+ +f+++ a elP +    l+G p ed e +f+ + ++  +  p  +  
      YP_700417.1  88 ELDPPEHTYYRQALQPLFSPKRMKELEPRIRDVINELIDDFAARGEAEFISEFAHELPTRVFLALMGWPLEDAE-MFTTTTDVALQGVPGGTee 180
                      ************************************************************************86.5666665554443333300 PP

  Cyp125(GramPos) 193 ..eedaaaasaellay.amklaeerkknPadd.ivtklveadi...dgek.lsddefgffvillavaGnettrnaithGmlafldnpdqWelyk 278
                        ++  +aa+ ++++y +  +a  r+ + + d + +++++a i   dg + l+d+e+  +  ll +aG  t + +++  ++ +++np q +   
      YP_700417.1 181 esATAREAAANQIFGYfGAIVAGVRSGEITSDtLTAQIINAPIemgDGVRlLTDEELYRMFFLLLIAGLHTVQGSLAWAIIHLANNPGQRQEIV 274
                      00334556677888883456677888777666155689999884344544278888887888899**********************9987655 PP

  Cyp125(GramPos) 279 r..erpetaadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarlei 370
                         +  ++a++ei+r  + v++  r a++d+e+ggv+i++g++++++  san d   f+dP+++ + r+pn hl+fG  G+h ciG++lar+e+
      YP_700417.1 275 DdpDTVSAAVEEILRIEAAVIAG-RRATRDVEIGGVTIREGDQLIVLLCSANRDGAEFEDPDELRVDRSPNRHLSFG-AGPHRCIGSHLARIEL 366
                      4124567899**********986.6699************************************************8.69************** PP

  Cyp125(GramPos) 371 dlifnaiadklPdlkkl.geperlrsgwlngikel 404
                      +l +++i ++lPd + + ++p  l+++ + g  +l
      YP_700417.1 367 KLAMEEIHKRLPDYRLVpEDPPILHATQVRGCIRL 401
                      *************9988356777777777666555 PP

>> YP_702911.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  181.1   0.0   3.4e-56   1.1e-53      97     395 ..      85     384 ..      11     397 .. 0.86

  Alignments for each domain:
  == domain 1  score: 181.1 bits;  conditional E-value: 3.4e-56
  Cyp125(GramPos)  97 lvllnkdapehtrlrkivsrlftPravealreeleerarkivkaa.aekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddp 189
                      ++ll+ d+p+h   r ++s+++ Pra+++lr++  + a+++v++  a++ ++d v+++a+ +Pl+   + +G+p++ re+l+ + ++ +++  p
      YP_702911.1  85 SILLEADPPHHDAPRAVLSKILGPRALQKLRAAWIQDAEALVDQLlANTTEFDAVTDLAAAFPLRVFPDAVGIPDAGRENLLPYGDHAFNAFGP 178
                      79****************************************9761667899**************************************9999 PP

  Cyp125(GramPos) 190 efaeed.aaaasaellayamklaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrerp 282
                           + +a   ael ++      +r +  +d   +++ +a   g +++ ++  + v  l  aG +tt n ++  + af+++pdqW   +++r 
      YP_702911.1 179 ANGLVEkGAPRVAELSGWVN-AQCARDALTGDGFGAQIWAAADRG-DITYEQAPLVVRSLLTAGVDTTVNGLAAVLYAFATHPDQWARLRENRT 270
                      98855525666777777654.344566666666666665555555.59******************************************9995 PP

  Cyp125(GramPos) 283 ..etaadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlif 374
                        +ta de vrw +Pv++f rta++dte+gg +i +g+++ +f  +an d+  +++Pe fd+ r+p  h+g+G  G h c+G+++arle + ++
      YP_702911.1 271 laRTAFDEAVRWESPVQTFFRTATRDTEIGGATIPDGKKILMFLGAANRDPRRWENPEVFDLGRNPSGHVGYG-MGIHQCVGQHVARLESEALL 363
                      4389********************************************************************7.7******************* PP

  Cyp125(GramPos) 375 naiadklPdlkklgeperlrs 395
                       a+a ++  l+ +g  +r   
      YP_702911.1 364 TALASRVHSLEIAGPVHRHLN 384
                      *********999998888655 PP

>> YP_702567.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  178.2   0.0   2.6e-55   8.1e-53      27     387 ..      26     379 ..      21     403 .. 0.82

  Alignments for each domain:
  == domain 1  score: 178.2 bits;  conditional E-value: 2.6e-55
  Cyp125(GramPos)  27 laelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddi.....ereqielqklvllnkdapehtrlrkivs 115
                      +++lr+t Pv  +  p++         W+vt ++ v+++  +   fss ++  i++++ ++     ++e   +   v++ +d+p+htr+r+ ++
      YP_702567.1  26 ITRLRETRPVSPMVFPDGHE------GWLVTGYDAVRQLMAD-TRFSSRQDIGILHVPYETpgmpaATEPSPQMPGVFIAMDPPDHTRMRRKLA 112
                      56677777776666655543......3999999999998765.5699999988888766541111123444455689***************** PP

  Cyp125(GramPos) 116 rlftPravealreeleerark.ivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgeddpefaeedaaaasaellayam 208
                       +ft + +++l++++ + a++ + + a+ +   d+v++ a  +P   i+ellGvp +dr++    s +++ +d+p    +d++aa   l +y  
      YP_702567.1 113 GAFTVKRMKQLEDHIIDVAERqLDAMARLTPPVDLVKEFALPVPSLVICELLGVPYADRDNFQVNSAKFLIKDQPL---DDKMAAYGALSTYLA 203
                      **************998887615567888899*****************************999999999998886...589999999****** PP

  Cyp125(GramPos) 209 klaeerkknPaddivtklveadidgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrer..petaadeivrwatPvtsfq 300
                      +l  ++++ P+ddi+++l + d    +l+ +e++    ll +aG ett n ++ G  a+l+np+q    + +      a++e++r+ +    f 
      YP_702567.1 204 DLVTRKRAAPGDDILSDLARDD----DLTIEELTGAAFLLLLAGHETTANMLALGAFALLENPEQLTELRTDPdlLPDAVEELLRYLSVADIFY 293
                      *****************96533....36677777677777889**********************98887764114579*************** PP

  Cyp125(GramPos) 301 rtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarleidlifnaiadklPdlkkl 387
                      r a+ed+elgg +i+ g  vv+   +an d++ fd+P+t+di r+   hl++G +G h c+G++larle++  f+ +  ++P l  +
      YP_702567.1 294 RYATEDIELGGETIRAGSTVVVSLLAANRDPQRFDNPDTLDIRRKARGHLSLG-HGVHLCLGQQLARLEMRAGFEGLLRRFPTLGLA 379
                      ***************************************************96.8*********************99999988543 PP

>> YP_700371.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  173.9   0.0   5.2e-54   1.6e-51       5     389 ..      17     393 ..      13     412 .. 0.85

  Alignments for each domain:
  == domain 1  score: 173.9 bits;  conditional E-value: 5.2e-54
  Cyp125(GramPos)   5 alpegfdftdpdllae.rlPveelaelrktaPvwWneqpeevggfddeGyWvvtkhkdvkevsrrsdvfsseentaivrfkddiereqielqkl 97 
                      +l  +fd+ dp l a   +  ++ a lr   P+ + ++         +G+W+vt+++d+  + r++++fss+ n  +             + k+
      YP_700371.1  17 DLMVDFDVYDPALAAPvDVFQAKAAGLRARGPILYSPH--------YGGHWIVTRYDDIFRILRDAETFSSYPNNLVD----------AGQGKF 92 
                      56678999999998752678888999999999999765........4799*********************9998773..........456899 PP

  Cyp125(GramPos)  98 vllnkdapehtrlrkivsrlftPravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelv....ged 187
                      + ++ d+peht+ r+ ++ lf P  +++l+ +++e  +++++  a++g+ +fv + a  lP +    l+G p +d e+  +w++  +    g+ 
      YP_700371.1  93 IPVEIDPPEHTQYRQALQPLFGPARMKELEPKIREIINELIDGFASSGRCEFVAEFAHALPTRVFLTLMGWPLDDAERFTEWTDIALqgipGAG 186
                      99*********************************************************************************86441111344 PP

  Cyp125(GramPos) 188 dpefaeedaaaasaellayamklaeerk..knPaddivtklveadi..dg..eklsddefgffvillavaGnettrnaithGmlafldnpdqWe 275
                      ++e  ++  a+a+++++ y  ++ ++ +  ++ ++ + ++++++ i  dg  + l+d+e+  +  ll +aG  t + ++   m+ +++np q +
      YP_700371.1 187 EEES-ARARAEAATNVFEYFGAFVDRVRsgQEDGESVTAQIINTPIlmDGveRYLTDEELRRMFFLLLIAGLHTVQGSLGWAMVHLANNPAQRN 279
                      4444.6667788899999987777654412455666777787775433552267*******99999**************************99 PP

  Cyp125(GramPos) 276 lykrerp..etaadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlar 367
                         ++      a++ei+r  + v s+ r a++d+e++gvk++ g+++ l+  san d+  fddP++f i r  n hl+fG  G+h c+G++lar
      YP_700371.1 280 ALIEDTSlvPDAVEEILRIEAAV-SMGRRATRDVEIAGVKVEAGDQLLLLLCSANRDDTEFDDPDAFVIDRGSNRHLSFG-AGPHRCLGSHLAR 371
                      877765412579*******9987.68999**************************************************8.69*********** PP

  Cyp125(GramPos) 368 leidlifnaiadklPdlkklge 389
                      le+ l +++i  ++Pd + +++
      YP_700371.1 372 LELTLAMEEIHRRIPDYALVES 393
                      ****************998864 PP

>> YP_702345.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  102.3   0.0   2.9e-32   9.1e-30     105     386 ..      82     378 ..      72     389 .. 0.86

  Alignments for each domain:
  == domain 1  score: 102.3 bits;  conditional E-value: 2.9e-32
  Cyp125(GramPos) 105 pehtrlrkivsrlftPravealreelee.rarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdw....snelvge...ddpe 190
                        h  lr +v  ++ P  v++  + l +  ar+ v++ +++g+ d+v+q +  + ++a+ +llG+ +   + l dw    sn+++++    d e
      YP_702345.1  82 EIHRDLRSMVDPALQPSEVDRWVDGLVRpIARRYVEQFENDGKADLVSQYCEPVSVRALGDLLGLNEVSSDTLRDWfhrlSNSFTNAgvdADGE 175
                      5699*****************9999865279**********************************99999988888555678877643226888 PP

  Cyp125(GramPos) 191 faeed....aaaasaellayamklaeerkknPaddivtklveadidgeklsddefgffvillav..aGnettrnaithGmlafldnpdqWelyk 278
                      f++ +    + +a+ae+ a    l ++ + +P d  +++ ++  +   ++ d e+ +  +++ +  a +e     ++  ++ +++ p+q e   
      YP_702345.1 176 FTNPEgfvqGDEAKAEIRAVVDPLIDKWTVHPDDSAISHWLHDGMPEGQVRDREYIYPTLFVYLlgAMQEP-GHGMASTLVGLFTRPEQLEAVI 268
                      88655333367899***************************99999889*****98855444440044554.467899999**********998 PP

  Cyp125(GramPos) 279 rerp..etaadeivrwatPvts.fqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyciGanlarle 369
                       e +    a+ e +rw++P+ s   r +++d+ lg+v + +g  v l y san+d  v+d P+++d++r p phl+fG +G h c G  +a+  
      YP_702345.1 269 DEPAliPRAISEGMRWTSPIWSaTARISTKDVTLGDVFLPEGSVVLLSYGSANHDTAVYDAPSDYDMTRPPLPHLAFG-SGNHACAGIYFANHV 361
                      87652246999*********75268****************************************************7.69************* PP

  Cyp125(GramPos) 370 idlifnaiadklPdlkk 386
                       ++ ++++ +++P+l++
      YP_702345.1 362 CRIGLEELFEAIPNLER 378
                      *************9986 PP

>> YP_703834.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   91.4   0.0     6e-29   1.9e-26      97     386 ..      75     379 ..      15     392 .. 0.83

  Alignments for each domain:
  == domain 1  score: 91.4 bits;  conditional E-value: 6e-29
  Cyp125(GramPos)  97 lvllnkdapehtrlrkivsrlftPravealreel.eerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvge... 186
                        +l  +  +h ++r+ +  ++ P av++  e+l + +ar++ +  ++ g+ d  e   + + ++++ +l+G+ +   e l +w ++l ++   
      YP_703834.1  75 GTVLAANGEQHEKIREWIDPQLRPSAVDSYVEALvRPQARSLLEGIEDLGAADIQEAYFAPISVRSVGDLMGLTEIPSETLVRWFETLAQSygn 168
                      567778889******************99777762679***********************************999999999999988753111 PP

  Cyp125(GramPos) 187 ....ddpefaeedaaaa....saellayamklaeerkknPaddivtklveadidgeklsd.de.fgffvillavaGnettrnaithGmlafldn 270
                          ++ +fa+  + +a    +ae++a    + ++ +++P + ++++ ++  + + ++ d +e +    ++l  a +e     +t  ++ ++++
      YP_703834.1 169 aevdENGNFANPGPFEAgdrvKAEIVAAVGPMLDHWTEHPDHTLISHWLHDGMPDGQVRDrSEiYPNIYVFLLGALQEPG-HVMTTTLAGLFQH 261
                      111133445554444432222789999999*****************998877555555514423334445555666655.679999******* PP

  Cyp125(GramPos) 271 pdqWelykrerp..etaadeivrwatPvts.fqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdpnphlgfGgtGahyci 361
                      pdq e    +      a++e  rw +P+ s   + a +d+ +gg  +  g  v+l y san de v+++ e+++i r   phl+fG  G h c 
      YP_703834.1 262 PDQLERVIDDPTliPRAVNEGARWVAPIWSaAVKVAGRDVTIGGIDLPTGTPVMLAYGSANRDESVWENAEAYEIDRPIMPHLAFG-AGNHACA 354
                      *****987765412579***********7615799**************************************************8.59***** PP

  Cyp125(GramPos) 362 GanlarleidlifnaiadklPdlkk 386
                      G+ l    +++ ++a+ +++P++++
      YP_703834.1 355 GTYLGTAIVRIALEALFETIPNIEP 379
                      **********************986 PP

>> YP_702742.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   42.1   0.0   5.4e-14   1.7e-11     304     369 ..       2      67 ..       1      83 [. 0.92

  Alignments for each domain:
  == domain 1  score: 42.1 bits;  conditional E-value: 5.4e-14
  Cyp125(GramPos) 304 ledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilrdp.nphlgfGgtGahyciGanlarle 369
                      l+dt lgg ++  g  + l++ +an d+ vf++P+++ + r     hl+fG +G h+c+ + lar+e
      YP_702742.1   2 LADTTLGGAELPAGSHLLLLWGAANRDPAVFERPDEIVLDRPHiRSHLAFG-KGVHFCVASALARME 67 
                      68***********************************999964379****7.8*************8 PP

>> YP_705149.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?    1.3   0.0      0.13        42      99     132 ..     104     137 ..     100     144 .. 0.84
   2 !   30.4   0.0   1.9e-10   5.9e-08     214     376 ..     230     420 ..     219     425 .. 0.77

  Alignments for each domain:
  == domain 1  score: 1.3 bits;  conditional E-value: 0.13
  Cyp125(GramPos)  99 llnkdapehtrlrkivsrlftPravealreelee 132
                      l+  d  eh   r+i++++ft   +++  e+l+ 
      YP_705149.1 104 LMLLDGDEHLAHRRIMQQAFTRDRLSRYTEALHP 137
                      4557999**************9999998888875 PP

  == domain 2  score: 30.4 bits;  conditional E-value: 1.9e-10
  Cyp125(GramPos) 214 rkknPaddivtklveadi.dgeklsddefgffvillavaGnettrnaithGmlafldnpdqWelykrerpe....t..............aade 288
                      r++  ++d+ + l + +  +g+++sdd+   ++i+l +a  +t+  +++  m+ + ++pd  e  +r+ ++    t              ++ e
      YP_705149.1 230 RRAGAGEDLFSALCHIESeEGQRFSDDDVVNHMIFLLMAAHDTSTITLSTMMQYLGQHPDWQERCRRQSAAlgtsTptydqldeltdldlVMKE 323
                      7778899******99876379****************************9999999***99899988765322211212211111111113457 PP

  Cyp125(GramPos) 289 ivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilr......dpnph.lg...fGgtGahyciGanlarleidl 372
                       +r   Pv  + r a+edte+ g  i +g  + ++ + + + +e + dPe+fd  r      ++  h  +   f g G h c+G ++a  ei+ 
      YP_705149.1 324 CLRLVPPVPVVARRAVEDTEVLGHYIPRGTYMSVVVHFTHHMPEYWPDPERFDPERfaperrEDKVHrFAwepF-GGGVHKCLGMHFAGAEIKT 416
                      89*****************************9999999999999999999999765322222455553232225.569**************98 PP

  Cyp125(GramPos) 373 ifna 376
                      ++++
      YP_705149.1 417 VMHH 420
                      8765 PP

>> YP_703037.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?    2.7   0.0     0.051        16     100     130 ..     107     137 ..     103     152 .. 0.82
   2 !   11.5   0.0   0.00011     0.035     213     275 ..     231     293 ..     223     303 .. 0.83
   3 !   18.1   0.0   1.1e-06   0.00034     286     376 ..     323     422 ..     317     426 .. 0.81

  Alignments for each domain:
  == domain 1  score: 2.7 bits;  conditional E-value: 0.051
  Cyp125(GramPos) 100 lnkdapehtrlrkivsrlftPravealreel 130
                      +  d+ eh r r+i++++ft   +e+  + l
      YP_703037.1 107 MLLDSEEHLRHRRIMQQAFTRSRLENTVDVL 137
                      567999**************88777666655 PP

  == domain 2  score: 11.5 bits;  conditional E-value: 0.00011
  Cyp125(GramPos) 213 erkknPaddivtklveadid.geklsddefgffvillavaGnettrnaithGmlafldnpdqWe 275
                      +r+++ +dd+ + l + + + g+++sdde   ++i+l +a  +t+  +i+  m+ + ++  qW+
      YP_703037.1 231 ARRSHETDDLFSVLCHIESEaGQRFSDDEVVDHMIFLLMAAHDTSTITISTMMQYLGQH-PQWQ 293
                      5888999*******999875499**********************99998766555555.5886 PP

  == domain 3  score: 18.1 bits;  conditional E-value: 1.1e-06
  Cyp125(GramPos) 286 adeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilr......dpnphl....gfGgtGahyciGanlarle 369
                      + e +r  +Pv  + r a++dt++ g  +  g   +++ + + + ++ + +Pe+fd  r      ++  h      fG  G h c+G ++a  e
      YP_703037.1 323 MKECLRIVSPVPVMARRAVRDTQVQGHFVPAGTYAAVVPHFTHHMPQYWPEPERFDPERfaehrrEDKVHRyawePFG-GGVHKCLGMHFAGAE 415
                      6799*******************************9999999999999999999986542222126666642222465.59************* PP

  Cyp125(GramPos) 370 idlifna 376
                      i+ i+++
      YP_703037.1 416 IKAIMHH 422
                      *998765 PP

>> YP_704615.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !    6.1   0.0    0.0045       1.4     198     273 ..     196     278 ..     185     287 .. 0.76
   2 !   26.7   0.0   2.5e-09   7.9e-07     283     377 ..     308     411 ..     300     431 .. 0.85

  Alignments for each domain:
  == domain 1  score: 6.1 bits;  conditional E-value: 0.0045
  Cyp125(GramPos) 198 aasaellayamklaeerkknP.....addivtklvea.didge.klsddefgffvillavaGnettrnaithGmlafldnpdq 273
                      +a ++l+a   ++ + r +nP       d++  lv+  d dg+ +++ de++ + i ++ aG  tt  + +  ++ +l +pd 
      YP_704615.1 196 EARVQLVALVQEIMNGRIENPpqgkeDRDMLDVLVSIkDEDGNeRFTADEITGMFISMMFAGHHTTSGTAAWTLIELLRHPDY 278
                      56667777777777777777633322358999998752556652799999988888999*************9*999999985 PP

  == domain 2  score: 26.7 bits;  conditional E-value: 2.5e-09
  Cyp125(GramPos) 283 etaadeivrwatPvtsfqrtaledtelggvkikkgqrvvlfyrsanfdeevfddPetfdilr..dpnph........lgfGgtGahyciGanla 366
                      e++  e +r   P++ + r a  + e+gg +i +++ v+   + +n   e f +P+tfd  r  dpn          + fG  G h c+Ga +a
      YP_704615.1 308 EAVLKETLRLHPPLIILLRVARGEFEVGGYRIAENDLVAATPAISNRIAEDFPNPDTFDPERyiDPNQEdivnrwtwIPFG-AGRHRCVGAAFA 400
                      56678999***************************************************9996677753332222225675.79********** PP

  Cyp125(GramPos) 367 rleidlifnai 377
                       ++++ if  +
      YP_704615.1 401 LMQLKAIFSIL 411
                      ********865 PP

>> YP_704571.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   -3.4   0.0       3.5   1.1e+03     104     148 ..     136     180 ..     133     196 .. 0.75
   2 !   14.5   0.0   1.3e-05    0.0041     213     377 ..     270     464 ..     242     475 .. 0.68

  Alignments for each domain:
  == domain 1  score: -3.4 bits;  conditional E-value: 3.5
  Cyp125(GramPos) 104 apehtrlrkivsrlftPravealreeleerarkivkaaaekgsgd 148
                      +p+ ++  ++++ +ft  a+++ + ++ + a +++++  e+  g 
      YP_704571.1 136 EPNWSKAHNLLAPAFTKSAMRSYHRTMLDVAGELTEHWDERVDGS 180
                      688889999999999999999999988888888877766665555 PP

  == domain 2  score: 14.5 bits;  conditional E-value: 1.3e-05
  Cyp125(GramPos) 213 erkknPaddivtklveadi..dgeklsddefgffvillavaGnettrnaithGmlafldnpd.....q......W.......elykrer.peta 285
                      + ++   +d++  +++a    d +++ + ++ + v+ + vaG ett  a++  +  +  +pd     q      W       e   + r  + +
      YP_704571.1 270 DSAEDGPEDLLELMLRAARenDPHRIDELNIRHQVVTFLVAGHETTSGALSFALYYLSRHPDvlakaQaevdavWgdeepafEQIAKLRyVRRV 363
                      445556678888888775511667889999999*********************9999999844444322222222221111222222213456 PP

  Cyp125(GramPos) 286 adeivrwatPvtsfqrtaledtelggv.kikkgqrvvlfyrsanfdeevfddPetfd.......ilrdpnphlg.fGgtGahyciGanlarlei 370
                       de +r      ++ r a+ dt l g+  +k g+ v ++  +   d+   ddPe+fd        +r+   h+    gtG   ciG ++a  e 
      YP_704571.1 364 LDESLRLWPTAPAYGREATVDTTLVGKyPMKVGDWVLVLIPALHRDPVWGDDPEAFDpdhflpeRIRSRPAHVYkPFGTGERACIGRQFALHES 457
                      788888777778999******99988726888999988888888887777889998744443333577777874234799999**999998887 PP

  Cyp125(GramPos) 371 dlifnai 377
                       l++  i
      YP_704571.1 458 VLVLGTI 464
                      7776655 PP

>> YP_703045.1  no_gene_name-cytochrome P450
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   13.3   0.0   3.1e-05    0.0096     142     188 ..      12      58 ..       1      77 [. 0.81

  Alignments for each domain:
  == domain 1  score: 13.3 bits;  conditional E-value: 3.1e-05
  Cyp125(GramPos) 142 aekgsgdfveqvavelPlqaiaellGvpqedreklfdwsnelvgedd 188
                      +++ ++d+v  +a+ +Pl+   + +G+p+  re+l+++ ++l+++  
      YP_703045.1  12 ETATTFDMVPALAAAFPLRVFPDAVGIPEVGRENLLSYGDHLFNAFG 58 
                      346789********************************999987655 PP

>> YP_702741.1  no_gene_name-hypothetical protein
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 ?   10.6   0.0    0.0002     0.064     121     172 ..       7      58 ..       2      60 .] 0.91

  Alignments for each domain:
  == domain 1  score: 10.6 bits;  conditional E-value: 0.0002
  Cyp125(GramPos) 121 ravealreeleerarkivkaaaekgsgdfveqvavelPlqaiaellGvpqed 172
                      r ++al  ++++   ++ ++  ++g+ d+ + +a +lP   +a+l+G+ ++d
      YP_702741.1   7 RRIRALTPTVHSLVDTLWDEGLRDGQIDWASAMANRLPPAMVAHLIGLTESD 58 
                      889999999999999999******************************9987 PP



Internal pipeline statistics summary:
-------------------------------------
Query model(s):                              1  (408 nodes)
Target sequences:                         7211  (2366103 residues searched)
Passed MSV filter:                       199  (0.0275967); expected 144.2 (0.02)
Passed bias filter:                      186  (0.0257939); expected 144.2 (0.02)
Passed Vit filter:                        30  (0.00416031); expected 7.2 (0.001)
Passed Fwd filter:                        23  (0.00318957); expected 0.1 (1e-05)
Initial search space (Z):               7211  [actual number of targets]
Domain search space  (domZ):              23  [number of targets reported over threshold]
# CPU time: 0.29u 0.02s 00:00:00.31 Elapsed: 00:00:00.11
# Mc/sec: 8776.09
//

'''