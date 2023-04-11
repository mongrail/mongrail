# finding the unique haplotypes from the set of haplotypes 
# observed in the colllection of both the populations as a whole

{
    array[$1]++
}

END{
    for(item in array)
    {
	print item
    }
}
