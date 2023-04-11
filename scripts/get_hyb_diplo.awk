# finding the occurrence of a specific haplotype
# from an already transposed individual data file
BEGIN {OFS = "\t"}

NR == 1 {next}

{
    # query_hap="11111"
    for(j = 2; j <= NF; j++)
    {
	if(j == 2)
	{
	    split($j,hap_init,"[|/]") # genotype 0|1 (marker 1)
	    hap_1=hap_init[1]      # hap_1 = 0
	    hap_2=hap_init[2]      # hap_2 = 1
	}
	else
	{
	    split($j,hap,"[|/]")       # genotype 1|1 (marker 2)
	    hap_1=hap_1""hap[1]     #hap_1 = 01
	    hap_2=hap_2""hap[2]     #hap_2 = 11
	}
    }

}

END{
    print hap_1, hap_2
}
