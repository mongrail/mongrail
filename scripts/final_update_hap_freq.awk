# using the haplotype counts to calculate the posterior mean (which is an estimate of the haplotype frequencies)

NR==FNR{
    uniq[$1]++;
    next
}

{
    if(uniq[$1]==1)
    {
	A[$1]=($2+(1.0/length(uniq)))/(1.0+$3)
	uniq[$1]=0
    }
    denominator=(1.0+$3)
}

END{
    for(item in uniq)
    {
	if(uniq[item]==1)
	{
	    A[item]=(1.0/length(uniq))/denominator
	}
    }
    for(hap in A)
    {
	printf ("%s:%f\t",hap, A[hap])
    }
    printf("\n")
}
