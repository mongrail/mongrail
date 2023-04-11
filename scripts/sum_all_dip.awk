# Sum over all possible diplotypes under a fixed model to obtain the log-probability of multilocus genotype of a hybrid individual under the six different models

NR != 1{
    for (i = 1; i <= 6; i++){
	value=exp($i)
	sum[i]= sum[i] + value
    }

}

END{
    for (item in sum){
	printf("%f\t",log(sum[item]))
    }
    printf("\n")
}
