# find the posterior probability for each scaffold from lklhd files and calculate the combined posterior
# probability at the end

BEGIN{
    printf("%s\t\t%-s\t%-s\t%-s\t%-s\t%-s\t%-s\n","Region","PostPr(a)","PostPr(b)","PostPr(c)","PostPr(d)","PostPr(e)","PostPr(f)")
    # printf("----------------------------------------\n\n")
}


{
    printf("%s\t",$1)
    denom=0
    for(i = 2; i <= NF; i++){
	denom=denom+exp($i)
	sum[i]=sum[i]+($i)
    }
    SUM=0
    for(i = 2; i <= NF; i++){
	p[i]=exp($i)/denom
	if(i==NF)
	{
	    probf=(1-SUM)
	    printf("%f\t",probf)
	}
	else
	{
	    SUM=SUM+p[i]
	    printf("%f\t",p[i])
	}
    }
    printf("\n")
}

END{
    # printf("-----------------------------------------\n\n")
    printf("%s\t\t","PostProb(all)")
    for (item in sum){
	combined_denom=combined_denom+exp(sum[item])
    }
    for (item in sum){
	combined_post_prob=exp(sum[item])/combined_denom
	printf("%f\t",combined_post_prob)
    }
    printf("\n")
}
