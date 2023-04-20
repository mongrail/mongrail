#!/bin/bash

function vcf_to_GT
{
    # local input_popA
    # local input_popB
    # local input_hyb
    local input_popA input_popB input_hyb scaffold_info
    local sample_popA sample_popB sample_hyb scaffold_info
    local header_popA header_popB header_hyb 
    
    input_popA="$1"
    input_popB="$2"
    input_hyb="$3"
    scaffold_info="$4"

    
    sample_popA=`bcftools query -l ${input_popA} | awk 'BEGIN{ORS="\t"} {print $1}'`
    sample_popB=`bcftools query -l ${input_popB} | awk 'BEGIN{ORS="\t"} {print $1}'`
    sample_hyb=`bcftools query -l ${input_hyb} | awk 'BEGIN{ORS="\t"} {print $1}'`
    # scaffold_info=`bcftools query -f '%CHROM\n' ${input_popA} | uniq -c`
    
    header_popA=`echo chrom:pos$'\t'${sample_popA}`                              
    header_popB=`echo chrom:pos$'\t'${sample_popB}`                              
    header_hyb=`echo chrom:pos$'\t'${sample_hyb}`
    
    while read -r n_markers scaffold_name                                        
    do                                                                           
	# echo "${n_markers} ****** ------- ${scaffold_name}"                      
	output_GT_popA="A_${scaffold_name}.GT"                                   
	output_GT_popB="B_${scaffold_name}.GT"                                   
	output_GT_hyb="H_${scaffold_name}.GT"                                    
        
	echo "${header_popA}" > ${output_GT_popA}                                
	echo "${header_popB}" > ${output_GT_popB}                                
	echo "${header_hyb}" > ${output_GT_hyb}                                  
        
	if [[ ${n_markers} -le 10 ]]                                             
	then                                                                     
	    bcftools query -r ${scaffold_name} -f '%CHROM:%POS\t\t[\t%GT]\n' ${input_popA} >> ${output_GT_popA}
	    bcftools query -r ${scaffold_name} -f '%CHROM:%POS\t\t[\t%GT]\n' ${input_popB} >> ${output_GT_popB}
	    bcftools query -r ${scaffold_name} -f '%CHROM:%POS\t\t[\t%GT]\n' ${input_hyb} >> ${output_GT_hyb}
	else                                                                     
	    echo "WARNING: More than 10 markers specified for ${scaffold_name}. Using the first 10 markers only!\n\n"
	    bcftools query -r ${scaffold_name} -f '%CHROM:%POS\t\t[\t%GT]\n' ${input_popA} >> ${output_GT_popA}
	    bcftools query -r ${scaffold_name} -f '%CHROM:%POS\t\t[\t%GT]\n' ${input_popB} >> ${output_GT_popB}
	    bcftools query -r ${scaffold_name} -f '%CHROM:%POS\t\t[\t%GT]\n' ${input_hyb} >> ${output_GT_hyb}
	    sed -i -n '1,11p' ${output_GT_popA}
	    sed -i -n '1,11p' ${output_GT_popB}
	    sed -i -n '1,11p' ${output_GT_hyb} 
	fi
    done < <(echo "${scaffold_info}")

}

function create_chrom
{
    local input_popA recom_flag 
    input_popA="$1"
    recom_flag="$2"
    
    if [[ ${recom_flag} == 0 ]]
    then
    	local input_recom_file 
    	input_recom_file="$3"
    	while read -r scaffold_name recom_rate                                         
    	do                                                                             
    	    markers=`bcftools query -r ${scaffold_name} -f '%POS\t' ${input_popA}`     
    	    n_marker=`echo "${markers}" | awk '{print NF}'`                            
    	    chrom_filename="${scaffold_name}.chrom"                                    
    	    if [[ ${n_marker} -le 10 ]]                                                
    	    then                                                                       
    		last_marker=`echo "${markers}" | awk '{print ($NF+(10**6))/(10**6)}'`  
    		chrom_length=`echo "${last_marker}/1" | bc`                            
    		chrom_info=`echo ${scaffold_name}$'\t'${chrom_length}$'\t'${recom_rate}$'\t'${n_marker}$'\t'${markers}`
    		echo "${chrom_info}" > ${chrom_filename}
    	    else                                                                       
    		truncated_n_marker=10
    		ten_markers=`echo "${markers}" | awk 'BEGIN{ORS="\t"}{for(i=1; i<=10; i++) print $i}'`
    		last_marker=`echo "${ten_markers}" | awk '{print ($NF+(10**6))/(10**6)}'`
    		chrom_length=`echo "${last_marker}/1" | bc`
    		chrom_info=`echo ${scaffold_name}$'\t'${chrom_length}$'\t'${recom_rate}$'\t'${truncated_n_marker}$'\t'${ten_markers}`
    		echo "${chrom_info}" > ${chrom_filename}                                                  
    	    fi
    	    sed -i "s/${scaffold_name}/1/g" ${chrom_filename}
    	done < <(echo "${input_recom_file}")
	
    else
    	local recom_rate scaffold_info
    	recom_rate="$3"
    	scaffold_info="$4"
    	# scaffold_info=`bcftools query -f '%CHROM\n' ${input_popA} | uniq`
    	while read -r scaffold_name                                        
    	do                                                                           
    	    markers=`bcftools query -r ${scaffold_name} -f '%POS\t' ${input_popA}`     
    	    n_marker=`echo "${markers}" | awk '{print NF}'`                            
    	    chrom_filename="${scaffold_name}.chrom"                                    
    	    if [[ ${n_marker} -le 10 ]]                                                
    	    then                                                                       
    		last_marker=`echo "${markers}" | awk '{print ($NF+(10**6))/(10**6)}'`  
    		chrom_length=`echo "${last_marker}/1" | bc`                            
    		chrom_info=`echo ${scaffold_name}$'\t'${chrom_length}$'\t'${recom_rate}$'\t'${n_marker}$'\t'${markers}`
    		echo "${chrom_info}" > ${chrom_filename}
    	    else                                                                       
    		truncated_n_marker=10
    		ten_markers=`echo "${markers}" | awk 'BEGIN{ORS="\t"}{for(i=1; i<=10; i++) print $i}'`
    		last_marker=`echo "${ten_markers}" | awk '{print ($NF+(10**6))/(10**6)}'`
    		chrom_length=`echo "${last_marker}/1" | bc`
    		chrom_info=`echo ${scaffold_name}$'\t'${chrom_length}$'\t'${recom_rate}$'\t'${truncated_n_marker}$'\t'${ten_markers}`
    		echo "${chrom_info}" > ${chrom_filename}                                                  
    	    fi
    	    sed -i "s/${scaffold_name}/1/g" ${chrom_filename}
    	done < <(echo "${scaffold_info}")
    fi
}

function create_pop_sim
{
    local scaffold_info
    scaffold_info="$1"
    while read -r scaffold_name
    do
	popA_scaffold_filename="A_${scaffold_name}.GT"                          
	popB_scaffold_filename="B_${scaffold_name}.GT"                          
	hybrid_scaffold_filename="H_${scaffold_name}.GT"


	header=`awk 'BEGIN {ORS="\t"}; {print $1}' ${hybrid_scaffold_filename}`        
	header+=$'\n'                                                                  
	number_hybrids=`awk '{print NF; exit}' ${hybrid_scaffold_filename}`            
	count=0                                                                        
	for (( j=2; j<=${number_hybrids}; j++))                                        
	do                                                                             
	    indv_header=`awk 'BEGIN {ORS="\t"}; {print $1}' ${hybrid_scaffold_filename}`
	    indv_header+=$'\n'                                                         
	    indv_count=0                                                               
            
	    hap1=`awk -v i=$j '{print $1, $i}' ${hybrid_scaffold_filename} | ./scripts/transpose.sh | awk -f ./scripts/get_hyb_diplo.awk | awk '{print $1}'`
	    hap2=`awk -v i=$j '{print $1, $i}' ${hybrid_scaffold_filename} | ./scripts/transpose.sh | awk -f ./scripts/get_hyb_diplo.awk | awk '{print $2}'`
	    hap_1=`echo "ibase=2;${hap1}" | bc`                                        
	    hap_2=`echo "ibase=2;${hap2}" | bc`                                        
            
	    if [ ${hap_1} != ${hap_2} ]                                                
	    then                                                                       
		all_dip=`./gendiplo 0 ${hap_1} ${hap_2}`                               
	    else                                                                       
		# reversed to be in parity with the gendiplo output format             
		rev_hap_1=`echo "${hap1}" | rev`                                       
		rev_hap_2=`echo "${hap2}" | rev`                                       
		all_dip=`echo "${rev_hap_1}" "${rev_hap_2}"`                           
	    fi                                                                         
            
	    while read line                                                            
	    do                                                                         
		count=$((count+1))                                                     
		indv_count=$((indv_count+1))                                           
		print_string=""                                                        
		haplo1=`echo "${line}" | awk '{print $1}'`                             
		haplo2=`echo "${line}" | awk '{print $2}'`                             
		hap1=`echo "${haplo1}" | rev`                                          
		hap2=`echo "${haplo2}" | rev`                                          
		# echo "i${count}: HAP1:${hap1} HAP2:${hap2}"                          
		for (( i=0; i<${#hap1}; i++))                                          
		do                                                                     
		    string="${hap1:$i:1}"                                              
		    string="${string}|"${hap2:$i:1}""                                  
		    print_string+=`echo ${string}$'\t'`                                
                    
		done                                                                   
		# echo i${count}$'\t'"${print_string}"                                 
		header+=`echo i${count}$'\t'"${print_string}"`                         
		header+=$'\n'                                                          
		
		indv_header+=`echo i${indv_count}$'\t'"${print_string}"`               
		indv_header+=$'\n'                                                     
		
	    done < <(echo "${all_dip}")                                                
            
	    indv_filename=`awk -v i=$j 'NR==1{print $i}' ${hybrid_scaffold_filename}`  
	    indv_filename="${indv_filename}_${scaffold_name}.sim"                      
	    indv_all_dip=`echo "${indv_header}" | ./scripts/transpose.sh`                         
	    echo "${indv_all_dip}" > "${indv_filename}"                                
	    sed -i "s/${scaffold_name}/1/g" ${indv_filename}                           
	done

	final_formatted_indv=`echo "${header}" | ./scripts/transpose.sh`
	
	
	
	
	# ${final_formatted_indv} contains the multi locus genotype data for hybrids (all possible phases compatible with the sampled diplotypes)
	# obtained after running the script "adjust_print_all_dip" on multi genotype data of hybrid individuals
	
	
	hyb=`echo "${final_formatted_indv}" | ./scripts/transpose.sh | awk -f ./scripts/hap_frequency.awk`  
	A_pop=`./scripts/transpose.sh ${popA_scaffold_filename} | awk -f ./scripts/hap_frequency.awk`       
	B_pop=`./scripts/transpose.sh ${popB_scaffold_filename} | awk -f ./scripts/hap_frequency.awk`       
	all_hap=`echo "${A_pop}" | awk '{print $1}'`                                   
	all_hap="${all_hap}"$'\n'"`echo "${B_pop}" | awk '{print $1}'`"                
	all_hap="${all_hap}"$'\n'"`echo "${hyb}" | awk '{print $1}'`"                  
	uniq_hap=`echo "${all_hap}" | awk -f ./scripts/uniq_hap.awk`                             
	A_pop_freq=`echo "1"$'\t'`                                                     
	B_pop_freq=`echo "1"$'\t'`                                                     
	A_pop_freq+=`awk -f ./scripts/final_update_hap_freq.awk <(echo "${uniq_hap}") <(echo "${A_pop}")`
	B_pop_freq+=`awk -f ./scripts/final_update_hap_freq.awk <(echo "${uniq_hap}") <(echo "${B_pop}")`
	
	A_pop_freq_filename="A_${scaffold_name}.popA"                                  
	B_pop_freq_filename="B_${scaffold_name}.popB"                                  
	
	echo "${A_pop_freq}" > ${A_pop_freq_filename}                                  
	echo "${B_pop_freq}" > ${B_pop_freq_filename}    

    done < <(echo "${scaffold_info}")
}


function create_lklhd
{
    local scaffold_info
    scaffold_info="$1"
    while read -r scaffold_name
    do
	A_freq_scaffold_filename="A_${scaffold_name}.popA"                      
	B_freq_scaffold_filename="B_${scaffold_name}.popB"
	hybrid_scaffold_filename="H_${scaffold_name}.GT"
	chrom_filename="${scaffold_name}.chrom"
	
	number_hybrids=`awk '{print NF; exit}' ${hybrid_scaffold_filename}`
	count=0
	for (( j=2; j<=${number_hybrids}; j++))
	do
	    
	    hybrid_name=`awk -v i=$j 'NR==1{print $i}' ${hybrid_scaffold_filename}`
	    hybrid_scaffold_sim_filename="${hybrid_name}_${scaffold_name}.sim"
	    hybrid_scaffold_out_filename="${hybrid_name}_${scaffold_name}.out"


	    ./mongrail -c ${chrom_filename}  -A ${A_freq_scaffold_filename} -B ${B_freq_scaffold_filename} -i ${hybrid_scaffold_sim_filename} -o ${hybrid_scaffold_out_filename} &&
		hybrid_lklhd_filename="${hybrid_name}.lklhd"
	    columns_A_to_F=`awk -f ./scripts/sum_all_dip.awk ${hybrid_scaffold_out_filename}`
	    output=`echo ${scaffold_name}$'\t'"${columns_A_to_F}"`
	    echo "${output}" >> ${hybrid_lklhd_filename}

	done
    done < <(echo "${scaffold_info}")

}


function create_PostProb
{
    local input_hyb
    input_hyb="$1"
    sample_hyb=`bcftools query -l ${input_hyb} | awk '{print $1}'`
    verbose_output="post_prob.txt"
    output="output.txt"
    post_prob_all_strings="Indiv PostPr(a) PostPr(b) PostPr(c) PostPr(d) PostPr(e) PostPr(f)"
    post_prob_all_strings+=$'\n'
    while read -r hybrid_name
    do
	hybrid_lklhd_filename="${hybrid_name}.lklhd"
	header="Indiv: ${hybrid_name}"
	echo -e "${header}\n" >> ${verbose_output}
	post_prob_output=`awk -f ./scripts/post_prob.awk ${hybrid_lklhd_filename}`
	echo "${post_prob_output}" | column -t >> ${verbose_output}
	echo -e "\n\n" >> ${verbose_output}
	last_line=`echo "${post_prob_output}" | tail -1 | sed "s/PostProb(all)/${hybrid_name}/g"`
	# echo "${last_line}" | column -t >> ${output}
	post_prob_all_strings+=`echo "${last_line}"`
	post_prob_all_strings+=$'\n'
    done < <(echo "${sample_hyb}")

    echo "${post_prob_all_strings}" | column -t > ${output}

}


while getopts ':A:B:i:r:R:h' opt; do                                                         
    case "$opt" in                                                                           
	A)                                                                                   
            popA="$OPTARG"                                                                   
            ;;                                                                               
        
	B)                                                                                   
            popB="$OPTARG"                                                                   
            ;;                                                                               
        
	i)                                                                                   
            hybrid="$OPTARG"                                                                 
            ;;

	r)
	    r_rate="$OPTARG"
	    rFlag=1
	    ;;

	R)
	    recom_file="$OPTARG"
	    ;;
        
	h)
            echo -e "mongrail v1.0.0 mongrail\nhttps://github.com/mongrail/mongrail\n"
            echo "Usage: ./$(basename $0) [-A FILENAME] [-B FILENAME] [-i FILENAME] [-r VALUE | -R FILENAME]"
            exit 0                                                                           
            ;;                                                                               
        
	:)
	    echo -e "mongrail v1.0.0 mongrail\nhttps://github.com/mongrail/mongrail\n"
            echo -e "ERROR: Option requires an argument.\nUsage: ./$(basename $0) [-A FILENAME] [-B FILENAME] [-i FILENAME] [-r VALUE | -R FILENAME]"
            exit 1                                                                           
            ;;                                                                               
        
	?)
	echo -e "mongrail v1.0.0 mongrail\nhttps://github.com/mongrail/mongrail\n"
	echo -e "ERROR: Invalid command option.\nUsage: ./$(basename $0) [-A FILENAME] [-B FILENAME] [-i FILENAME] [-r VALUE | -R FILENAME]"
	exit 1                                                                               
	;;                                                                                   
    esac                                                                                     
done
# echo $OPTIND                                                                             

if ((OPTIND < 9))                                                                          
then                                                                                       
    
    if ((OPTIND == 1))                                                                     
    then
        echo -e "mongrail v1.0.0 mongrail\nhttps://github.com/mongrail/mongrail\n"
        echo -e "ERROR: No options specified!\nUsage: ./$(basename $0) [-A FILENAME] [-B FILENAME] [-i FILENAME] [-r VALUE | -R FILENAME]"
    else
        echo -e "mongrail v1.0.0 mongrail\nhttps://github.com/mongrail/mongrail\n"
        echo -e "ERROR: All options not specified!\nUsage: ./$(basename $0) [-A FILENAME] [-B FILENAME] [-i FILENAME] [-r VALUE | -R FILENAME]"
    fi                                                                                     
    exit 1                                                                                 
fi                       


no_marker_scaffold=`bcftools query -f '%CHROM\n' ${popA} | uniq -c`
# scaffold_name=`echo "${no_marker_scaffold}" | awk 'BEGIN{ORS="\n"}{print $2}'`
chrom_name=`bcftools query -f '%CHROM\n' ${popA} | uniq`
# echo "${scaffold_name}"


vcf_to_GT ${popA} ${popB} ${hybrid} "${no_marker_scaffold}"    
create_chrom ${popA} ${rFlag} ${r_rate} "${chrom_name}"
create_pop_sim "${chrom_name}"
create_lklhd "${chrom_name}"
create_PostProb ${hybrid}


mkdir tmp
mv *.chrom ./tmp/
mv *.popA ./tmp/
mv *.popB ./tmp/
mv *.GT ./tmp/
mv *.sim ./tmp/
mv *.out ./tmp/
mv *.lklhd ./tmp/

echo ""
