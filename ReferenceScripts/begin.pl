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
				print "\n-------------------------------------------------------------------\n";
				@d=split(/\s+/,$c[$i]);
				print $q."\t".$c[0]."\t$d[6]\t$d[7]\t$d[8]\t$d[10]\t$d[11]\n" if $d[6]<1;
				print "\n";
				print $q;
				print "\n";
				print $c[0];
				print "\n";
				print  "\t#\t\tscore\tbias\tc-Evalue\ti-Evalue\thmmfrom\thmmto\t\talifrom\talito\t\tenvfrom\tenvto\t\tacc\n";
				print join("  \t",@d)."\n";
				print  "\t1\t2\t3\t4\t5\t\t6\t\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\n";
				print 'd[6] = '."$d[6]\n";
				print 'd[7] = '."$d[7]\n";
				print 'd[8] = '."$d[8]\n";
				print 'd[10] = '."$d[10]\n";
				print 'd[11] = '."$d[11]\n";
			}
		}
		@a=();
	}
	else
	{
		push(@a,$_);
	}
}

