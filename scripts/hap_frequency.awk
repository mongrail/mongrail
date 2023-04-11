# find unique haplotypes for a single population and their corresponding counts
# from a transposed individual data file

# BEGIN {OFS = "\t"}

NR == 1 {next}

{
    for(j = 2; j <= NF; j++)
    {
	if(j == 2)
	{
	    split($j,hap_init,"|") # genotype 0|1 (marker 1)
	    hap_1=hap_init[1]      # hap_1 = 0
	    hap_2=hap_init[2]      # hap_2 = 1
	}
	else
	{
	    split($j,hap,"|")       # genotype 1|1 (marker 2)
	    hap_1=hap_1""hap[1]     #hap_1 = 01
	    hap_2=hap_2""hap[2]     #hap_2 = 11
	}
    }
    array[hap_1]++
    array[hap_2]++
}

END{
    for(item in array)
    {
	# printf ("%s:%f\t", item, (array[item]/((NR-1)*2)))
	# print item, array[item], ((NR-1)*2), (array[item]/((NR-1)*2))
	printf ("%s\t%d\t%d\t%f\n", item, array[item], ((NR-1)*2), (array[item]/((NR-1)*2)))
    }
    # printf ("\n")
}
