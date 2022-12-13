#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<glib.h>
#include<regex.h>
#include<limits.h>
#include<stddef.h>
#include<unistd.h>
#include<ctype.h>
#include<stdbool.h>
#include<assert.h>

#define MAX_CHROMOSOMES 23
#define TOTAL_MARKERS 500
#define BUFFER_SIZE 60000
#define MAX_HAPSEQS 1024
#define MAX_INDIVIDUALS 11000
#define MAX_FILENAME 100
#define MAX_INDIVIDUAL_NAME 10
#define NO_MODELS 6
#undef DEBUG

GError *error = NULL;
GIOChannel *infile_chrom = NULL;
GIOChannel *infile_A = NULL;
GIOChannel *infile_B = NULL;
GIOChannel *infile_indv = NULL;
GIOChannel *outfile = NULL;
GIOChannel *outfile_AT = NULL;

typedef struct chromosome_all_info
{
  unsigned int chrom_id;
  double chrom_length;
  double chrom_recom_rate;
  unsigned int n_loci;
  double *markers;
  GHashTable* hashA;
  GHashTable* hashB;
  unsigned int *haplotype_1;
  unsigned int *haplotype_2;
  struct chromosome_all_info *next;
} chrom_data;


chrom_data* make_data_node(unsigned int index, double length, double recom_rate, unsigned int no_loci, double *loci)
{
  chrom_data *new_data_node;
  if((new_data_node = malloc(sizeof(chrom_data))) == NULL)
    { fprintf(stderr, "Oops, out of memory!"); exit(1);}

  new_data_node->chrom_id = index;
  new_data_node->chrom_length = length;
  new_data_node->chrom_recom_rate = recom_rate;
  new_data_node->n_loci = no_loci;
  new_data_node->markers = malloc(no_loci * sizeof(double));
  assert(new_data_node->markers != NULL);
  for(int i = 0; i < no_loci; i++)
    {
      new_data_node->markers[i] = loci[i];
    }
  new_data_node->hashA = NULL;
  new_data_node->hashB = NULL;
  new_data_node->haplotype_1 = NULL;
  new_data_node->haplotype_2 = NULL;
  new_data_node->next = NULL;
  return(new_data_node);
}

void string_to_markers(char *string, unsigned int no_loci, double **markers)
{
  char *second_regexString = "([0-9]+)";
  char *second_string;                                                
  second_string = string; 
  regex_t *second_regexCompiled = (regex_t*)malloc(sizeof(regex_t));
  assert(second_regexCompiled != NULL);
  regmatch_t *pmatch;                                                 
  if(regcomp(second_regexCompiled, second_regexString, REG_EXTENDED|REG_NEWLINE))             
    {                                                                 
      printf("Could not compile regular expression.\n");              
      exit(1);                                                        
    }                                                                 
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*second_regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0))  
    {                                                                 
      printf("Nothing matched with ""%s""\n", second_string);         
      exit(1);                                                        
    }   
  unsigned int n_loci = 0;                                            
  double *loci_position;                                              
  loci_position = malloc(TOTAL_MARKERS * sizeof(double));
  assert(loci_position != NULL);
  do {                                                                
    if (pmatch[0].rm_so != -1) {        /* The regex is matching part\
					   of a string */                                                       
      char *submatch;                                                 
      double val;                                                     
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;            
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch, second_string + pmatch[0].rm_so, matchlen+1); 
      submatch[matchlen]='\0';                                        
      val = atof(submatch);                                           
      loci_position[n_loci] = val;                                    
      free(submatch);                                                 
    };                                                                
    second_string += pmatch[0].rm_eo;   /* Restart from last match */ 
    n_loci++;                                                         
  } while(!regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch,0));
  *markers = malloc(n_loci * sizeof(double));
  assert(markers != NULL);
  for(int i = 0; i < n_loci; i++ )                                    
    {
      *(*markers + i) = loci_position[i];
#ifdef DEBUG
      printf("Locus[%d]: %lf\t %lf\n", i+1, loci_position[i],*(*markers + i));
#endif
    }  
  regfree(second_regexCompiled);                                      
  free(second_regexCompiled);                                         
  free(pmatch);
  free(loci_position);

}

void read_chrom(char *filename, GIOChannel *infile, unsigned int **chromosome_index, unsigned int **chromosome_no_loci, double **chromosome_length, double **chromosome_recom_rate, unsigned int *no_chromosomes, double **all_markers, unsigned int *no_total_markers)
{
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char * regexString = "([0-9]+)[[:blank:]]+([0-9]+\\.?[0-9]*)[[:blank:]]+([0-9\
]+\\.?[0-9]*)[[:blank:]]+([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 6;
                                                                                
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  unsigned int *chrom_index, *chrom_no_loci;
  double *chrom_length, *chrom_recom_rate;
  double *chrom_markers;
  chrom_index = calloc(MAX_CHROMOSOMES, sizeof(unsigned int));
  chrom_no_loci = calloc(MAX_CHROMOSOMES, sizeof(unsigned int));
  chrom_length = calloc(MAX_CHROMOSOMES, sizeof(double));
  chrom_recom_rate = calloc(MAX_CHROMOSOMES, sizeof(double));
  chrom_markers = calloc(MAX_CHROMOSOMES * TOTAL_MARKERS, sizeof(double));
  assert(chrom_index != NULL);
  assert(chrom_no_loci != NULL);
  assert(chrom_length != NULL);
  assert(chrom_recom_rate != NULL);
  assert(chrom_markers != NULL);


  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };

  unsigned int c = 0;
  unsigned int chrom_no = 0;
  unsigned int val_1, val_4;
  double val_2, val_3;
  int marker_index = 0;
  for(c = 0; c < MAX_CHROMOSOMES; c++)
    {
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	{
	  unsigned int g = 0;
	  unsigned int offset = 0;
	  char sourceCopy[strlen(source) + 1];
	  strcpy(sourceCopy, source);
	  sourceCopy[groupArray[g].rm_eo] = 0;
	  for (g = 0; g < maxGroups; g++)
	    {
	      if (groupArray[g].rm_so == (size_t)-1)
		break;  // No more groups

	      if(g == 0)
		{
		  offset = groupArray[g].rm_eo;
		}
	      if(g == 1)
		{
		  char *submatch_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
		  val_1 = atoi(submatch_1);
		  chrom_index[chrom_no] = val_1;
#ifdef DEBUG
		  printf("Chromosome ID: %u\n\n", val_1);
#endif
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  char *submatch_2;
		  size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_2 = (char*)malloc(matchlen_2+1);
		  assert(submatch_2 != NULL);
		  strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		  submatch_2[matchlen_2]='\0';
		  val_2 = atof(submatch_2);
		  chrom_length[chrom_no] = val_2;
#ifdef DEBUG
		  printf("Chromosome length (in Mb): %lf\n\n", val_2);
#endif
		  free(submatch_2);
		}
	      if(g == 3)
		{
		  char *submatch_3;
		  size_t matchlen_3 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_3 = (char*)malloc(matchlen_3+1);
		  assert(submatch_3 != NULL);
		  strncpy(submatch_3, sourceCopy + groupArray[g].rm_so, matchlen_3+1);
		  submatch_3[matchlen_3]='\0';
		  val_3 = atof(submatch_3);
		  chrom_recom_rate[chrom_no] = val_3;
#ifdef DEBUG
		  printf("Recombination rate (in cM/Mb): %lf\n\n", val_3);
#endif
		  free(submatch_3);
		}
	      if(g == 4)
		{
		  char *submatch_4;
		  size_t matchlen_4 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_4 = (char*)malloc(matchlen_4+1);
		  assert(submatch_4 != NULL);
		  strncpy(submatch_4, sourceCopy + groupArray[g].rm_so, matchlen_4+1);
		  submatch_4[matchlen_4]='\0';
		  val_4 = atoi(submatch_4);
		  chrom_no_loci[chrom_no] = val_4;
#ifdef DEBUG
		  printf("Number of markers: %u\n\n", val_4);
#endif
		  free(submatch_4);
		}
	      if(g == 5)
	      	{
		  char *submatch_5;
		  size_t matchlen_5 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_5 = (char*)malloc(matchlen_5+1);
		  assert(submatch_5 != NULL);
		  strncpy(submatch_5, sourceCopy + groupArray[g].rm_so, matchlen_5+1);
		  submatch_5[matchlen_5]='\0';
		  /* printf("%s,\n\n",submatch_5); */
		  double *loci;
		  string_to_markers(submatch_5, chrom_no_loci[chrom_no], &loci);
		  for(int i = 0; i < val_4; i++)
		    {
		      chrom_markers[marker_index + i] = loci[i];
		    }
		  free(loci);
		  free(submatch_5);
	      	}
	    }
	  marker_index = marker_index + val_4;	  
	  source = source + offset;
	  chrom_no++;
	}
    }
  /* printf("# chrom:%u\n\n\n", chrom_no); */
  
  *chromosome_index = malloc(chrom_no * sizeof(unsigned int));
  *chromosome_no_loci = malloc(chrom_no * sizeof(unsigned int));
  *chromosome_length = malloc(chrom_no * sizeof(double));
  *chromosome_recom_rate = malloc(chrom_no * sizeof(double));
  *all_markers = malloc(marker_index * sizeof(double));
  assert(chromosome_index != NULL);
  assert(chromosome_no_loci != NULL);
  assert(chromosome_length != NULL);
  assert(chromosome_recom_rate != NULL);
  assert(all_markers != NULL);

  for(int i = 0; i < chrom_no; i++)
    {
      *(*chromosome_index + i) = chrom_index[i];
      *(*chromosome_no_loci + i) = chrom_no_loci[i];
      *(*chromosome_length + i) = chrom_length[i];
      *(*chromosome_recom_rate + i) = chrom_recom_rate[i];
    }
  for(int i = 0; i < marker_index; i++)
    {
      *(*all_markers + i) = chrom_markers[i];
    }

  *no_chromosomes = chrom_no;
  *no_total_markers = marker_index;
  regfree(&regexCompiled);
  g_free(complete_file);
  free(chrom_index);
  free(chrom_no_loci);
  free(chrom_length);
  free(chrom_recom_rate);
  free(chrom_markers);
}

chrom_data* data_linked_list_creation(unsigned int *chromosome_index, unsigned int *chromosome_no_loci, double *chromosome_length, double *chromosome_recom_rate, unsigned int no_chromosomes, double *chromosome_markers)
{
  chrom_data* head, *current;
  int first = 0;
  int marker_index = 0;
  for(int i = 0; i < no_chromosomes; i++)
    {
      double *loci;
      loci = chromosome_markers + marker_index;
      chrom_data *next_data_node = make_data_node(chromosome_index[i], chromosome_length[i], chromosome_recom_rate[i], chromosome_no_loci[i], loci);
      marker_index = marker_index + chromosome_no_loci[i];
      if(first == 0)
	{
	  head = next_data_node;
	  current = head;
	  first = 1;
	}
      else
	{
	  current->next = next_data_node;
	  current = current->next;
	}
    }
  return(head);

}


void hap_info_extract(char* string, unsigned int **sequence, double **frequency, unsigned int *no_sequence)
{
  char *second_string;
  second_string = string;
  char *second_regexString = "([0-1]+:[0]*\\.[0-9]+)";
  regex_t *second_regexCompiled = (regex_t*)malloc(sizeof(regex_t));
  assert(second_regexCompiled != NULL);
  regmatch_t *pmatch;
  if(regcomp(second_regexCompiled, second_regexString, REG_EXTENDED|REG_NEWLINE))
    {                                                                            
      printf("Could not compile regular expression.\n");
      exit(1);
    }                                                                            
  pmatch = (regmatch_t*)malloc(sizeof(regmatch_t)*second_regexCompiled->re_nsub);
  assert(pmatch != NULL);
  if(regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0))
    {                                                                                       
      printf("Nothing matched with ""%s""\n", second_string);
      exit(1);
    }
  unsigned int n_hapseq = 0;
  unsigned int *hap_seq;
  double *hap_freq;
  hap_seq = malloc(MAX_HAPSEQS * sizeof(unsigned int));
  hap_freq = malloc(MAX_HAPSEQS * sizeof(double));
  assert(hap_seq != NULL);
  assert(hap_freq != NULL);

  do {                                                                          
    if (pmatch[0].rm_so != -1) {        /* The regex is matching part of a string */
      char *submatch;
      size_t matchlen = pmatch[0].rm_eo - pmatch[0].rm_so;
      submatch = (char*)malloc(matchlen+1);
      assert(submatch != NULL);
      strncpy(submatch, second_string + pmatch[0].rm_so, matchlen+1);
      submatch[matchlen]='\0';
      /* printf("%s\n",submatch); */
      char *source;
      source = submatch;
      char * regexString = "([0-1]+):([0]*\\.[0-9]+)";
      size_t maxGroups = 3;
      regex_t regexCompiled;
      regmatch_t groupArray[maxGroups];
      if (regcomp(&regexCompiled, regexString, REG_EXTENDED))                                      
        {                                                                                          
          printf("Could not compile regular expression.\n");
          exit(1);
        };
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
        {
          unsigned int g = 0;
          for (g = 0; g < maxGroups; g++)                                                          
            {                                                                                      
              if (groupArray[g].rm_so == (size_t)-1)                                               
                break;  // No more groups
              char sourceCopy[strlen(source) + 1];
              strcpy(sourceCopy, source);
              sourceCopy[groupArray[g].rm_eo] = 0;
	      if(g == 1)
		{
		  char *submatch_1;
		  long val_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
		  val_1 = strtol(submatch_1, NULL, 2);
		  hap_seq[n_hapseq] = (unsigned int) val_1;
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  char *submatch_2;
		  double val_2;
		  size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_2 = (char*)malloc(matchlen_2+1);
		  assert(submatch_2 != NULL);
		  strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		  submatch_2[matchlen_2]='\0';
		  val_2 = atof(submatch_2);
		  hap_freq[n_hapseq] = val_2;
		  free(submatch_2);
		}
            }
        }
      regfree(&regexCompiled);
      free(submatch);
    };
    second_string += pmatch[0].rm_eo;   /* Restart from last match */
    n_hapseq++;
  } while(!regexec(second_regexCompiled, second_string, second_regexCompiled->re_nsub, pmatch, 0)); 
  *sequence = malloc(n_hapseq * sizeof(unsigned int));
  *frequency = malloc(n_hapseq * sizeof(double));
  assert(sequence != NULL);
  assert(frequency != NULL);

  for(int i = 0; i < n_hapseq; i++)
    {
      *(*sequence + i)  = hap_seq[i];
      *(*frequency + i) = hap_freq[i];
    }
  *no_sequence = n_hapseq;
  regfree(second_regexCompiled);
  free(second_regexCompiled);
  free(pmatch);
  free(hap_seq);
  free(hap_freq);
}

void iterator(gpointer key, gpointer value, gpointer user_data)  {
  printf(user_data, *(guint32*) key, *(gdouble*)value);
}

void read_pop(char *filename, GIOChannel *infile, chrom_data *head, char pop_name)
{
  chrom_data* current;
  current = head;
  char *complete_file;
  if(g_io_channel_read_to_end(infile,&complete_file,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *source;
  source = complete_file;
  char * regexString = "([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 3;
                                                                                
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };

  unsigned int chrom_no = 0;
  unsigned int *temporary_key;
  gdouble* temporary_freq;
  while(current != NULL)
    {
      if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	{
	  unsigned int g = 0;
	  unsigned int offset = 0;
	  char sourceCopy[strlen(source) + 1];
	  strcpy(sourceCopy, source);
	  sourceCopy[groupArray[g].rm_eo] = 0;
	  for (g = 0; g < maxGroups; g++)
	    {
	      if (groupArray[g].rm_so == (size_t)-1)
		break;  // No more groups

	      if(g == 0)
		{
		  offset = groupArray[g].rm_eo;
		}
	      if(g == 1)
		{
		  char *submatch_1;
		  size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		  submatch_1 = (char*)malloc(matchlen_1+1);
		  assert(submatch_1 != NULL);
		  strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		  submatch_1[matchlen_1]='\0';
#ifdef DEBUG
		  printf("%s\n\n",submatch_1);
#endif
		  free(submatch_1);
		}
	      if(g == 2)
		{
		  if(pop_name == 'A')
		    {		 
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
#ifdef DEBUG
		      printf("~~ %s ~~\n\n", submatch_2);
#endif
		      unsigned int *hap_sequence, no_hap_sequence;
		      double *hap_frequency;
		      hap_info_extract(submatch_2, &hap_sequence, &hap_frequency, &no_hap_sequence);
#ifdef DEBUG
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  printf("Sequence: %u\tFrequency: %lf\n",hap_sequence[i], hap_frequency[i]);
			}
#endif
		      current->hashA = g_hash_table_new/* _full */(g_direct_hash, g_direct_equal/* , free, free */);
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  temporary_key = g_new(guint,1);
			  *temporary_key = hap_sequence[i];
			  temporary_freq = g_new(gdouble,1);
			  *temporary_freq = hap_frequency[i];
			  g_hash_table_insert(current->hashA, temporary_key, temporary_freq);
			}
		      
		      free(hap_sequence);
		      free(hap_frequency);
		      free(submatch_2);
		      
		    }
		  else
		    {
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
#ifdef DEBUG
		      printf("~~ %s ~~\n\n", submatch_2);
#endif
		      unsigned int *hap_sequence, no_hap_sequence;
		      double *hap_frequency;
		      hap_info_extract(submatch_2, &hap_sequence, &hap_frequency, &no_hap_sequence);
#ifdef DEBUG
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  printf("Sequence: %u\tFrequency: %lf\n",hap_sequence[i], hap_frequency[i]);
			}
#endif
		      current->hashB = g_hash_table_new/* _full */(g_direct_hash, g_direct_equal/* , free, free */);
		      for(int i = 0; i < no_hap_sequence; i++)
			{
			  temporary_key = g_new(guint,1);
			  *temporary_key = hap_sequence[i];
			  temporary_freq = g_new(gdouble,1);
			  *temporary_freq = hap_frequency[i];
			  g_hash_table_insert(current->hashB, temporary_key, temporary_freq);
			}
		      		      
		      free(hap_sequence);
		      free(hap_frequency);
		      free(submatch_2);

		    }

		}
	    }
	  source = source + offset;
	  chrom_no++;
	  /* printf("\n\n"); */
	}
      current = current->next;
    }

  regfree(&regexCompiled);
  g_free(complete_file);
}


void update_haplotype_indv(char *SNP, unsigned int *hap1, unsigned int *hap2, unsigned int position)
{
  char *second_regexString = "(0|1)\\|(0|1)";
  char *second_string;
  second_string = SNP;
#ifdef DEBUG
  printf("%s\n\n",second_string);
#endif
  size_t max_subgroup = 3;
  regex_t second_regexCompiled;
  regmatch_t pmatch[max_subgroup];
  unsigned int m = 0;
  char *cursor;
  
  if(regcomp(&second_regexCompiled, second_regexString, REG_EXTENDED))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    }
  cursor = second_string;
  /* unsigned int flag; */


  while(regexec(&second_regexCompiled, cursor, max_subgroup, pmatch, 0) == 0)
    {
#ifdef DEBUG
      flag = regexec(&second_regexCompiled, cursor, max_subgroup, pmatch, 0);
      printf("~~ Flag:%u ~~\n",flag);
#endif
      unsigned int g1 = 0;
      unsigned int offset_1 = 0;
      for (g1 = 0; g1 < max_subgroup; g1++)
	{
	  if (pmatch[g1].rm_so == (size_t)-1)
	    break;  // No more groups

	  if (g1 == 0)
	    offset_1 = pmatch[g1].rm_eo;

	  char cursorCopy[strlen(cursor) + 1];
	  strcpy(cursorCopy, cursor);
	  cursorCopy[pmatch[g1].rm_eo] = 0;
#ifdef DEBUG
	  printf("Match %u, Group %u: [%2u-%2u]: %s\n",
		 m, g1, pmatch[g1].rm_so, pmatch[g1].rm_eo,		 
		 cursorCopy + pmatch[g1].rm_so);
#endif
	  if(g1 == 1)
	  	{
	  	  unsigned int val_1;
	  	  val_1 = (unsigned int) atoi(cursorCopy + pmatch[g1].rm_so);
		  /* printf("%u\n",val_1); */
	  	  hap1[m] = (hap1[m] |(val_1 << position));
	  	}
	  if(g1 == 2)
	  	{
	  	  unsigned int val_2;
	  	  val_2 = (unsigned int) atoi(cursorCopy + pmatch[g1].rm_so);
		  /* printf("%u\n",val_2); */
	  	  hap2[m] = (hap2[m] |(val_2 << position));
	  	}

	}
      cursor += offset_1;
      m++;
    }

  regfree(&second_regexCompiled);

}


void update_linked_list_indv_data(char *individual_data, chrom_data *head, unsigned int no_indv)
{
  chrom_data *current;
  current = head;
  /* printf("~~~~\n\n"); */
  /* printf("%s\n", individual_data); */
  char *source;
  source = individual_data;
  char * regexString = "([0-9]+):([0-9]+)[[:blank:]]+(.+)";
  size_t maxGroups = 4;
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];

  if (regcomp(&regexCompiled, regexString, REG_EXTENDED|REG_NEWLINE))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    };
  
  while(current != NULL)
    {
      unsigned int *hap_1, *hap_2;
      hap_1 = calloc(no_indv, sizeof(unsigned int));
      hap_2 = calloc(no_indv, sizeof(unsigned int));
      assert(hap_1 != NULL);
      assert(hap_2 != NULL);
      for(int i = 0; i < current->n_loci; i++)
	{
	  if (regexec(&regexCompiled, source, maxGroups, groupArray, 0) == 0)
	    {
	      unsigned int g = 0;
	      unsigned int offset = 0;
	      char sourceCopy[strlen(source) + 1];
	      strcpy(sourceCopy, source);
	      sourceCopy[groupArray[g].rm_eo] = 0;

	      for (g = 0; g < maxGroups; g++)
		{
		  if (groupArray[g].rm_so == (size_t)-1)
		    break;  // No more groups

		  if(g == 0)
		    {
		      offset = groupArray[g].rm_eo;
		    }
		  if(g == 1)
		    {
		      char *submatch_1;
		      size_t matchlen_1 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_1 = (char*)malloc(matchlen_1+1);
		      assert(submatch_1 != NULL);
		      strncpy(submatch_1, sourceCopy + groupArray[g].rm_so, matchlen_1+1);
		      submatch_1[matchlen_1]='\0';
		      /* printf("%s\n\n",submatch_1); */
		      free(submatch_1);
		    }
		  if(g == 2)
		    {
		      char *submatch_2;
		      size_t matchlen_2 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_2 = (char*)malloc(matchlen_2+1);
		      assert(submatch_2 != NULL);
		      strncpy(submatch_2, sourceCopy + groupArray[g].rm_so, matchlen_2+1);
		      submatch_2[matchlen_2]='\0';
		      /* printf("~~ %s ~~\n\n", submatch_2); */
		      free(submatch_2);
		    }
		  if(g == 3)
		    {
		      char *submatch_3;
		      size_t matchlen_3 = groupArray[g].rm_eo - groupArray[g].rm_so;
		      submatch_3 = (char*)malloc(matchlen_3+1);
		      assert(submatch_3 != NULL);
		      strncpy(submatch_3, sourceCopy + groupArray[g].rm_so, matchlen_3+1);
		      submatch_3[matchlen_3]='\0';
		      /* printf("~~ %s ~~\n\n", submatch_3); */
		      update_haplotype_indv(submatch_3, hap_1, hap_2, current->n_loci-i-1);
		      free(submatch_3);
		    }
		}
	      source = source + offset;
	    }
	}
      current->haplotype_1 = malloc(no_indv * sizeof(unsigned int));
      current->haplotype_2 = malloc(no_indv * sizeof(unsigned int));
      assert(current->haplotype_1 != NULL);
      assert(current->haplotype_2 != NULL);
      for(int j = 0; j < no_indv; j++)
	{
	  current->haplotype_1[j] = hap_1[j];
	  current->haplotype_2[j] = hap_2[j];
#ifdef DEBUG
	  printf("%u\t%u\n",hap_1[j],hap_2[j]);
#endif
	}
      /* printf("\n\n"); */
      free(hap_1);
      free(hap_2);
      current = current->next;
    }
  regfree(&regexCompiled);

}

void read_multiple_indv(char *filename, GIOChannel *infile, chrom_data *head, unsigned int *n_indv/* , char *indv_names */)
{
  char *header;
  if(g_io_channel_read_line(infile,&header,NULL,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }

  char *regexString = "([[:alnum:]:_]+)";
  size_t max_subgroup = 1;
  regex_t regexCompiled;
  regmatch_t pmatch[max_subgroup];
  unsigned int m = 0;
  char *cursor;

  /* char *individual_ID[MAX_INDIVIDUALS]; */
  
  if(regcomp(&regexCompiled, regexString, REG_EXTENDED))
    {
      printf("Could not compile regular expression.\n");
      exit(1);
    }
  cursor = header;
  /* unsigned int flag; */
  while(regexec(&regexCompiled, cursor, max_subgroup, pmatch, 0) == 0)
    {
#ifdef DEBUG
      flag = regexec(&regexCompiled, cursor, max_subgroup, pmatch, 0);
      printf("~~ Flag:%u ~~\n",flag);
#endif
      unsigned int g1 = 0;
      unsigned int offset = 0;
      for (g1 = 0; g1 < max_subgroup; g1++)
        {
          if (pmatch[g1].rm_so == (size_t)-1)
            break;  // No more groups

          if (g1 == 0)
            offset = pmatch[g1].rm_eo;

          char cursorCopy[strlen(cursor) + 1];
          strcpy(cursorCopy, cursor);
          cursorCopy[pmatch[g1].rm_eo] = 0;
          /* printf("Match %u, Group %u: [%2u-%2u]: %s\n", */
          /*        m, g1, pmatch[g1].rm_so, pmatch[g1].rm_eo, */
          /*        cursorCopy + pmatch[g1].rm_so); */
	  /* if(m != 0){ */
	  /*   strcpy(individual_ID[m],cursorCopy + pmatch[g1].rm_so); */
	  /* } */
        }
      cursor += offset;
      m++;
    }

  /* printf("Count regexec: %u\n",m); */


  
  
  char *complete_file_without_header;
  if(g_io_channel_read_to_end(infile,&complete_file_without_header,NULL,&error) != G_IO_STATUS_NORMAL)
    {
      fprintf(stderr,"Found the file: '%s' but could not read the rest of the line\n ", filename);
      exit(1);
    }
  char *data;
  data = complete_file_without_header;
  /* printf("\n\n%s\n\n",data); */
  *n_indv = (m-1);

  update_linked_list_indv_data(data,head,m-1);

  regfree(&regexCompiled);
  g_free(complete_file_without_header);
  g_free(header);

}

double prob_change(double interval_length, double recom_rate)
{
  double lambda, value;
  
  lambda = (interval_length * recom_rate);
  value = 1 - (exp(-lambda) * cosh(lambda));
  return(value);
}

void all_marker_prob_change(double **inter_mark_prob_change,double **inter_marker_length, double recom_rate, unsigned int n_loci)
{
  for(int i = 0; i < n_loci; i++)
    {
      *(*inter_mark_prob_change + i) = prob_change(*(*inter_marker_length + i),recom_rate);
    }  
}


unsigned int returnBit(unsigned int n, unsigned int k)
{
  unsigned int temp;
  temp = n;
  temp = (temp & (1 << k));
  if(temp == 0)
    {
      return(0);
    }
  else
    {
      return(1);
    }
}

double Q_component(unsigned int start_state, unsigned int ancestry_state,double **inter_marker_prob_change, unsigned int n_loci, double pi)
{
  unsigned int previous_state, present_state;
  previous_state = start_state;
  double log_product;
  log_product = log(pi);
  for(int i = 0; i < n_loci; i++)
    {
      present_state = returnBit(ancestry_state, n_loci - i -1);
      if(present_state != previous_state)
	{
	  log_product = log_product + log(*(*inter_marker_prob_change + i));
	}
      else
	{
	  log_product = log_product + log(1 - *(*inter_marker_prob_change + i));
	}
      previous_state = present_state;
    }
  double Q;
  Q = log_product;
  return(Q);
}

double Q(unsigned int ancestry_state, double **inter_marker_prob_change, unsigned int n_loci, double pi_0, double pi_1)
{
  double log_Q_0, log_Q_1;
  log_Q_0 = Q_component(0, ancestry_state, inter_marker_prob_change, n_loci, pi_0);
  log_Q_1 = Q_component(1, ancestry_state, inter_marker_prob_change, n_loci, pi_1);
  double scale;
  scale = abs(fmax(log_Q_0,log_Q_1)) - 0.5;
  
  double log_Q;
  log_Q = log(exp(log_Q_0 + scale) + exp(log_Q_1 + scale)) - scale;
  return(log_Q);
  
}


typedef struct hap
{
  unsigned int interval_ancestry_state;
  unsigned int identity_state;
  unsigned int sub_hap;
  double sub_hap_freq;
  struct hap *next;
} hap_ancestry;

unsigned int setBit(unsigned int n, unsigned int k)
{
  unsigned int temp;
  temp = n;
  temp = (temp | (1 << k));
  return(temp);
}

hap_ancestry* makeNode(unsigned int interval_ancestry, unsigned int identity, unsigned int interval_hap, double interval_hap_freq)
{
  hap_ancestry* newNode;
  if((newNode = malloc(sizeof(hap_ancestry))) == NULL)
    { fprintf(stderr, "Oops, out of memory!"); exit(1); }
  
  newNode->interval_ancestry_state = interval_ancestry;
  newNode->identity_state = identity;
  newNode->sub_hap = interval_hap;
  newNode->sub_hap_freq = interval_hap_freq;
  newNode->next = NULL;
  return(newNode);
}

void update_SNP_Linked_List(hap_ancestry* head, unsigned int SNP_hap)
{
  hap_ancestry* current;
  current =  head;
  while(current != NULL)
    {
      current->sub_hap = (current->identity_state & SNP_hap);
      current = current->next;
    }
}

hap_ancestry* Linked_List_creation(unsigned int ancestry_state, unsigned int SNP_hap, unsigned int n_loci)
{
  hap_ancestry* head;
  hap_ancestry* current;
  unsigned int previous_state, present_state, match;
  unsigned int n = 0;
  unsigned int no_nodes = 0;
  if((returnBit(ancestry_state, n_loci-1)) == 0)
    {
      previous_state = 0;
      match = setBit(n, n_loci-1);
    }
  else
    {
      previous_state = 1;
      match = setBit(n, n_loci-1);
    }

  int first = 0;
  for(int i = 1; i < n_loci; i++)
    {
      present_state = returnBit(ancestry_state, n_loci-i-1);
      if(present_state != previous_state)
	{
	  hap_ancestry* nextNode = makeNode(previous_state, match, 0, 0);
	  no_nodes = no_nodes + 1;
	  if(first == 0)
	    {
	      head = nextNode;
	      current = head;
	      first = 1;
	    }
	  else
	    {
	      current->next = nextNode;
	      current = current->next;
	    }
	  match = setBit(0, n_loci-i-1);
	  previous_state = present_state;
	}
      match = setBit(match, n_loci-i-1);
    }
  if(no_nodes != 0)
    {
      hap_ancestry* lastNode = makeNode(previous_state, match, 0, 0);
      current->next = lastNode;
      current = current->next;
    }
  else
    {
      head = makeNode(previous_state, match, 0, 0);
      current = head;
      current->next = NULL;
    }
  update_SNP_Linked_List(head, SNP_hap);
  return(head);

}



void delete_linked_list(hap_ancestry *head)
{
  hap_ancestry *temp;
  while(head != NULL)
    {
      temp = head;
      head = head->next;
      free(temp);
    }
}

unsigned long int merge(unsigned int identity_state, unsigned int sub_haplotype, unsigned int no_loci)
{
  unsigned long int temp;
  temp = (((unsigned long int)sub_haplotype << no_loci) | (unsigned long int)identity_state);
  return(temp);
}

void iterator_sub(gpointer key, gpointer value, gpointer user_data) {
  printf(user_data, *(guint64*)key, *(gdouble*)value);
}

void prob_SNP_given_z(unsigned int identity_state, unsigned int SNP,
		      chrom_data* data_node, GHashTable* hashA_subSNP,
		      GHashTable* hashB_subSNP,double *log_prob) 
{
  chrom_data* data_node_copy;
  data_node_copy = data_node;  
  double log_freq = 0;
  /* printf("\n\nLinked List Creation\n\n"); */
  /* printf("z = %u\tSNP = %u\tn_loci = %u\n",identity_state, SNP, data_node_copy->n_loci); */

  hap_ancestry* head;
  head = Linked_List_creation(identity_state, SNP, data_node_copy->n_loci);
  hap_ancestry* current;

  /* current = head; */
  /* while(current != NULL) */
  /*   { */
  /*     printf("ancestry state: %u\tidentity state: %u\tsub haplotype: %u\t frequency: %lf\n", current->interval_ancestry_state, current->identity_state, current->sub_hap, current->sub_hap_freq); */
  /*     current = current->next; */
  /*   } */

  /* g_hash_table_foreach(data_node_copy->hashA, (GHFunc)iterator, "Haplotype Sequence: %u has Frequency %lf \n"); */
  /* printf("\n\n\n"); */
  /* g_hash_table_foreach(data_node_copy->hashB, (GHFunc)iterator, "Haplotype Sequence: %u has Frequency %lf \n"); */

  
  /* printf("Ancestry state:%u\n",identity_state); */
  /* printf("SNP:%u\n",SNP); */

  unsigned int table_size_A, table_size_B;
  table_size_A = g_hash_table_size(data_node->hashA);
  table_size_B = g_hash_table_size(data_node->hashB);

  gpointer *keys_A, *keys_B;
  keys_A = g_hash_table_get_keys_as_array(data_node->hashA, &table_size_A);
  keys_B = g_hash_table_get_keys_as_array(data_node->hashB, &table_size_B);
  /* for(int i = 0; i < table_size_A; i++) */
  /*   { */
  /*     printf("%d:%u\n",i,*(guint*)keys_A[i]); */
  /*     printf("Freq:%lf\n",*(double*)g_hash_table_lookup(data_node_copy->hashA,keys_A[i])); */

  /*   } */


  unsigned long int *temporary_key;
  gdouble* temporary_freq;
  current = head;
  while(current != NULL)
    {
      if(current->interval_ancestry_state == 0)
	{
	  temporary_key = malloc(sizeof(unsigned long int));
	  assert(temporary_key != NULL);
	  *temporary_key = merge(current->identity_state, current->sub_hap, data_node_copy->n_loci);
	  /* printf("A:temp_key: %lu\n",*temporary_key); */
	  if(g_hash_table_contains(hashA_subSNP, temporary_key) == FALSE)
	    {

	      for(int i = 0; i < table_size_A; i++)
		{
		  if((*(guint*)keys_A[i] & current->identity_state) == current->sub_hap)
		    {
		      current->sub_hap_freq = current->sub_hap_freq + *(double*)g_hash_table_lookup(data_node_copy->hashA,keys_A[i]);
		  
		    }
		}

	      temporary_freq = g_new(gdouble,1);
	      *temporary_freq = current->sub_hap_freq;
	      g_hash_table_insert(hashA_subSNP, temporary_key, temporary_freq);
	    }
	  else
	    {
	      /* printf("PREVIOUSLY CALCULATED frequency of subSNP: %u given interval_identity_state: %u for Population A \n\n", current->sub_hap, current->identity_state); */

	      current->sub_hap_freq = *(gdouble*)g_hash_table_lookup(hashA_subSNP, temporary_key);
	    }
	}
      else
	{
	  temporary_key = malloc(sizeof(unsigned long int));
	  assert(temporary_key != NULL);
	  *temporary_key = merge(current->identity_state, current->sub_hap, data_node_copy->n_loci);
	  /* printf("B:temp_key: %lu\n",*temporary_key); */
	  if(g_hash_table_contains(hashB_subSNP, temporary_key) == FALSE)
	    {
	      for(int i = 0; i < table_size_B; i++)
		{
		  if((*(guint*)keys_B[i] & current->identity_state) == current->sub_hap)
		    {
		      current->sub_hap_freq = current->sub_hap_freq + *(double*)g_hash_table_lookup(data_node_copy->hashB,keys_B[i]);
		    }
		}

	      temporary_freq = g_new(gdouble,1);
	      *temporary_freq = current->sub_hap_freq;
	      g_hash_table_insert(hashB_subSNP, temporary_key, temporary_freq);
	    }
	  else
	    {
	      /* printf("PREVIOUSLY CALCULATED frequency of subSNP: %u given interval_identity_state: %u for Population B \n\n", current->sub_hap, current->identity_state); */
	      
	      current->sub_hap_freq = *(gdouble*)g_hash_table_lookup(hashB_subSNP, temporary_key);
	    }
	}
      
      current = current->next;
    }

  /* printf("\n\nFinal Updated Linked List with frequency of the sub_SNP haplotype information\n");   */
  current = head;
  while(current != NULL)
    {
      log_freq = (log_freq + log(current->sub_hap_freq));
      /* printf("ancestry state: %u\tidentity state: %u\tsub haplotype: %u\t frequency: %lf\n", current->interval_ancestry_state, current->identity_state, current->sub_hap, current->sub_hap_freq);       */
      current = current->next;
    }

  /* printf("\n\nThe probability of SNP haplotype %u given ancestral state %u is (in log_e):\t%lf\n\n", SNP, identity_state, log_freq);   */
  *log_prob = log_freq;

  g_free(keys_A);
  g_free(keys_B);

  delete_linked_list(head);

}

double recom_hap_log_prob(unsigned int SNP, double *Q_array, unsigned int no_ancestry_state,
			  chrom_data* data_node, GHashTable* hashA_subSNP,
			  GHashTable* hashB_subSNP/* , unsigned int no_hapseqA, */
			  /* unsigned int no_hapseqB, unsigned int *hapseq_A, */
			  /* unsigned int *hapseq_B */)
{
  chrom_data* data_node_copy;
  data_node_copy = data_node;

  double scale, maximum, temporary;
  double *log_p, *log_R;
  log_p = malloc(no_ancestry_state * sizeof(double));
  log_R = malloc(no_ancestry_state * sizeof(double));
  assert(log_p != NULL);
  assert(log_R != NULL);
  double final_log_prob;
  double sum = 0;
  
  prob_SNP_given_z(0,SNP,data_node_copy,hashA_subSNP,hashB_subSNP,/* no_hapseqA,no_hapseqB, */
		   /* hapseq_A,hapseq_B, */ &log_p[0]);
  log_R[0] = log_p[0] + Q_array[0];
  maximum = log_R[0];
  temporary = log_R[0];
  for(unsigned int i = 1; i < no_ancestry_state; i++)
    {
      prob_SNP_given_z(i,SNP,data_node_copy,hashA_subSNP,hashB_subSNP,/* no_hapseqA,no_hapseqB, */
		       /* hapseq_A,hapseq_B, */ &log_p[i]);
      log_R[i] = log_p[i] + Q_array[i];

      temporary = log_R[i];
      if(maximum < temporary)
	{
	  maximum = temporary;
	}
    }
  /* printf("%lf\n",maximum); */
  scale = abs(maximum)-0.5;
  for(unsigned int i = 0; i < no_ancestry_state; i++)
    {
      sum = sum + exp(log_R[i] + scale) ;
    }
  final_log_prob = log(sum) - scale;
  /* printf("\n\nLog(base e) of the recombinant haplotype: %lf\n\n", final_log_prob); */
  free(log_p);
  free(log_R);
  return(final_log_prob);
}


void model_inference(int individual_index, double *Q_array, unsigned int no_ancestry_state,
		     chrom_data* data_node, GHashTable* hashA_subSNP, GHashTable* hashB_subSNP,
		     double *indv_log_prob)
{
  if((data_node->haplotype_1[individual_index]^data_node->haplotype_2[individual_index]) == 0) /* Bitwise XOR between haplotypes; 0 means they are identical*/
    {
      /* Likelihood calculation under model A when haplotype configurations are identical */
      double log_p1_a, log_p2_a;
      prob_SNP_given_z(no_ancestry_state-1,data_node->haplotype_1[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p1_a);  
      prob_SNP_given_z(no_ancestry_state-1,data_node->haplotype_2[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p2_a);
      indv_log_prob[0] = (log_p1_a + log_p2_a);

      /* Likelihood calculation under model D when haplotype configurations are identical */
      double log_p1_d, log_p2_d;
      prob_SNP_given_z(0,data_node->haplotype_1[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p1_d);
      prob_SNP_given_z(0,data_node->haplotype_2[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p2_d);
      indv_log_prob[3] = (log_p1_d + log_p2_d);

      /* Likelihood calculation under model C when haplotype configurations are identical */
      indv_log_prob[2] = (log_p1_a + log_p2_d);

      double recombinant_haplotype_log_probability_1, recombinant_haplotype_log_probability_2;
      recombinant_haplotype_log_probability_1 = recom_hap_log_prob(data_node->haplotype_1[individual_index],
								   Q_array,no_ancestry_state,data_node,hashA_subSNP,hashB_subSNP);
  
      recombinant_haplotype_log_probability_2 = recom_hap_log_prob(data_node->haplotype_2[individual_index],
								   Q_array,no_ancestry_state,data_node,hashA_subSNP,hashB_subSNP);

      /* Likelihood calculation under model E when haplotype configurations are identical */
      indv_log_prob[4] = (log_p1_a + recombinant_haplotype_log_probability_2);

      /* Likelihood calculation under model B when haplotype configurations are identical */
      indv_log_prob[1] = (log_p1_d + recombinant_haplotype_log_probability_2);

      /* Likelihood calculation under model F when haplotype configurations are identical */
      indv_log_prob[5] = (recombinant_haplotype_log_probability_1 + recombinant_haplotype_log_probability_2);
    }
  else
    {
      /* Likelihood calculation under model A */
      double log_p1_a, log_p2_a;
      prob_SNP_given_z(no_ancestry_state-1,data_node->haplotype_1[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p1_a);
      prob_SNP_given_z(no_ancestry_state-1,data_node->haplotype_2[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p2_a);
      indv_log_prob[0] = (log(2) + log_p1_a + log_p2_a);

      /* Likelihood calculation under model D */
      double log_p1_d, log_p2_d;
      prob_SNP_given_z(0,data_node->haplotype_1[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p1_d);
      prob_SNP_given_z(0,data_node->haplotype_2[individual_index],
		       data_node,hashA_subSNP,hashB_subSNP,&log_p2_d);
      indv_log_prob[3] = (log(2) + log_p1_d + log_p2_d);


      /* Likelihood calculation under model C */
      double scale_c, first_term_c, second_term_c;
      first_term_c = (log_p1_a + log_p2_d);
      second_term_c = (log_p1_d + log_p2_a);
      scale_c = abs(fmax(first_term_c, second_term_c))-0.5;
      indv_log_prob[2] = log(exp(first_term_c + scale_c) + exp(second_term_c + scale_c))-scale_c;

      double recombinant_haplotype_log_probability_1, recombinant_haplotype_log_probability_2;
      recombinant_haplotype_log_probability_1 = recom_hap_log_prob(data_node->haplotype_1[individual_index],
								   Q_array,no_ancestry_state,data_node,hashA_subSNP,hashB_subSNP);
  
      recombinant_haplotype_log_probability_2 = recom_hap_log_prob(data_node->haplotype_2[individual_index],
								   Q_array,no_ancestry_state,data_node,hashA_subSNP,hashB_subSNP);


      /* Likelihood calculation under model E */
      double scale_e, first_term_e, second_term_e;
      first_term_e = (log_p1_a + recombinant_haplotype_log_probability_2);
      second_term_e = (log_p2_a + recombinant_haplotype_log_probability_1);
      scale_e = abs(fmax(first_term_e, second_term_e))-0.5;
      indv_log_prob[4] = log(exp(first_term_e + scale_e) + exp(second_term_e + scale_e))-scale_e;


      /* Likelihood calculation under model B */
      double scale_b, first_term_b, second_term_b;
      first_term_b = (log_p1_d + recombinant_haplotype_log_probability_2);
      second_term_b = (log_p2_d + recombinant_haplotype_log_probability_1);
      scale_b = abs(fmax(first_term_b, second_term_b))-0.5;
      indv_log_prob[1] = log(exp(first_term_b + scale_b) + exp(second_term_b + scale_b))-scale_b;

      /* Likelihood calculation under model F */
      indv_log_prob[5] = (log(2) + recombinant_haplotype_log_probability_1 + recombinant_haplotype_log_probability_2);


    }

#ifdef DEBUG  
  printf("The Log Likelihood under model (a) is: %lf\n", indv_log_prob[0]);
  printf("The Log Likelihood under model (b) is: %lf\n", indv_log_prob[1]);
  printf("The Log Likelihood under model (c) is: %lf\n", indv_log_prob[2]);
  printf("The Log Likelihood under model (d) is: %lf\n", indv_log_prob[3]);
  printf("The Log Likelihood under model (e) is: %lf\n", indv_log_prob[4]);
  printf("The Log Likelihood under model (f) is: %lf\n", indv_log_prob[5]);
#endif

}






void update_model_indv_prob(chrom_data *data_node, double **matrix, unsigned int no_individual)
{
  /* chrom_data *data_node_copy; */
  /* data_node_copy = data_node; */


  double *loci_position;
  loci_position = malloc(data_node->n_loci * sizeof(double));
  assert(loci_position != NULL);
  for(int i = 0; i < data_node->n_loci; i++)
    {
      loci_position[i] = (double) (data_node->markers[i] / (data_node->chrom_length * pow(10,6)));
    }

  double *interval_length;
  interval_length = malloc(data_node->n_loci * sizeof(double));
  assert(interval_length != NULL);
  for(int i = 0; i < data_node->n_loci; i++)
    {
      if(i == 0)
  	{
  	  interval_length[i] = loci_position[i];
  	}
      else
	{
	  interval_length[i] = loci_position[i] - loci_position[i-1];
	}
    }

  double r; /* r is the expectated number of recombination for a given length of chromosome */
  r = (data_node->chrom_recom_rate * data_node->chrom_length / 100.0);
  
  /* printf("The probability of state change for each inter-marker interval is given by\n"); */
  double *P;
  P = malloc(data_node->n_loci * sizeof(double));
  assert(P != NULL);
  all_marker_prob_change(&P,&interval_length,r,data_node->n_loci);
#ifdef DEBUG
  for(int i = 0; i < data_node->n_loci; i++)
    {
      printf("P[%d] = %lf\n",i,P[i]);
    }
#endif
  /* printf("\n\n"); */
  /* printf("The value of Q for each ancestry state is as follows:\n"); */
  unsigned int index, mid_index;
  index = pow(2,data_node->n_loci);
  mid_index = pow(2,data_node->n_loci - 1);
  double *Q_ancestry;
  Q_ancestry = malloc(index * sizeof(double));
  assert(Q_ancestry != NULL);
  for(unsigned int i = 0; i < mid_index; i++)
    {
      Q_ancestry[i] = Q(i, &P, data_node->n_loci, 0.5, 0.5);
      Q_ancestry[index-i-1] = Q_ancestry[i];
    }
#ifdef DEBUG
  for(unsigned int i = 0; i < index; i++)
    {
      printf("Q[%u] = %lf\n",i,Q_ancestry[i]);
    }
#endif
  
  GHashTable* hash_A_sub = g_hash_table_new_full(g_direct_hash, g_direct_equal, free, free);
  GHashTable* hash_B_sub = g_hash_table_new_full(g_direct_hash, g_direct_equal, free, free);

  
  for(int i = 0; i < no_individual; i++)
    {
      model_inference(i,Q_ancestry,index,data_node,hash_A_sub,hash_B_sub,matrix[i]);

    }


  g_hash_table_destroy(hash_A_sub);
  g_hash_table_destroy(hash_B_sub);
  free(interval_length);
  free(P);
  free(Q_ancestry);
  free(loci_position);

}


void largest(double array[], int n, GHashTable *hash_table,
	     char *best_model, char *next_best_model, double *bayes_factor)
{
  double largest, second_largest;
  int max_index, second_max_index;
  
  // Initialize largest and second largest element

  if(array[0] > array[1])
    {
      largest = array[0];
      second_largest = array[1];
      max_index = 0;
      second_max_index = 1;
    }
  else
    {
      largest = array[1];
      second_largest = array[0];
      max_index = 1;
      second_max_index = 0;
    }
  
  // Traverse array elements from third and
  // compare every element with current max and update second largest value
  for (int i = 2; i < n; i++)
    {
      if (array[i] > largest)
	{
	  second_largest = largest;
	  second_max_index = max_index;
	  largest = array[i];
	  max_index = i;
	}
      else
	{
	  if (array[i] > second_largest)
	    {
	      second_largest = array[i];
	      second_max_index = i;
	    }
	}
    }

  *best_model = *(gchar*)g_hash_table_lookup(hash_table,&max_index);
  *next_best_model = *(gchar*)g_hash_table_lookup(hash_table,&second_max_index);

  *bayes_factor = (largest/second_largest);
}


void free_a_hash_table_entry(gpointer key, gpointer value, gpointer user_data)
{
  g_free(key);
  g_free(value);
}

void free_all_key_value_entries (GHashTable *table)
{
    g_hash_table_foreach (table, free_a_hash_table_entry, NULL);
    g_hash_table_destroy (table);
}

void delete_chrom_data(chrom_data *head)
{
  chrom_data *temp;
  while(head != NULL)
    {
      temp = head;
      head = head->next;
      free(temp->markers);
      free_all_key_value_entries (temp->hashA);
      free_all_key_value_entries (temp->hashB);
      free(temp->haplotype_1);
      free(temp->haplotype_2);
      /* g_hash_table_destroy(temp->hashA); */
      /* g_hash_table_destroy(temp->hashB); */
      free(temp);
    }
}


int main(int argc, char **argv)
{
  int c;
  int pflag = 0, Aflag = 0, Bflag = 0, iflag = 0;
  char *filenames[4];
  for(int i = 0; i < 4; i++)
    {
      filenames[i] = calloc(MAX_FILENAME, sizeof(char));
      assert(filenames[i] != NULL);
    }
  while((c = getopt(argc, argv, ":p:A:B:i:")) != -1)
    {
      switch(c)
	{
	case 'p':
	  strcpy(filenames[0],optarg);
	  pflag = 1;
	  break;
	case 'A':
	  strcpy(filenames[1],optarg);
	  Aflag = 1;
	  break;
	case 'B':
	  strcpy(filenames[2],optarg);
	  Bflag = 1;
	  break;
	case 'i':
	  strcpy(filenames[3],optarg);
	  iflag = 1;
	  break;
	case ':':
	  if(optopt == 'p'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'A'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'B'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  else if(optopt == 'i'){
	    fprintf(stderr,"Missing argument for option -%c\n", optopt);}
	  for(int i = 0; i < 4; i++)
	    {
	      free(filenames[i]);
	    }
	  exit(1);
	case '?':
	  if (isprint (optopt))
	    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
	  else
	    fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
	  for(int i = 0; i < 4; i++)
	    {
	      free(filenames[i]);
	    }
	  exit(1);
	}
    }
  
  if(pflag == 0 || Aflag == 0 || Bflag == 0 || iflag == 0)
    {
      if(pflag == 0)
	{
	  fprintf(stderr,"%s: missing -p option\n",argv[0]);
	}
      if(Aflag == 0)
	{
	  fprintf(stderr,"%s: missing -A option\n",argv[0]);
	}
      if(Bflag == 0)
	{
	  fprintf(stderr,"%s: missing -B option\n",argv[0]);

	}
      if(iflag == 0)
      	{
      	  fprintf(stderr,"%s: missing -i option\n",argv[0]);
      	}
      for(int i = 0; i < 4; i++)
	{
	  free(filenames[i]);
	}
      exit(1);
    }

  
  if(pflag)
    {
      infile_chrom = g_io_channel_new_file(filenames[0],"r",&error);
      if(infile_chrom == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[0]);
	  exit(1);
	}
    }
  if(Aflag)
    {
      infile_A = g_io_channel_new_file(filenames[1],"r",&error);
      if(infile_A == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[1]);
	  exit(1);
	}

    }
  if(Bflag)
    {
      infile_B = g_io_channel_new_file(filenames[2],"r",&error);
      if(infile_B == NULL)
	{
	  fprintf(stderr,"Could not open file %s\n", filenames[2]);
	  exit(1);
	}
    }
  if(iflag)
    {
      infile_indv = g_io_channel_new_file(filenames[3],"r",&error);
      if(infile_indv == NULL)
  	{
  	  fprintf(stderr,"Could not open file %s\n", filenames[3]);
  	  exit(1);
  	}

    }


  unsigned int *chrom_index, *chrom_no_loci;
  unsigned int no_chrom;
  double *chrom_length, *chrom_recom_rate;
  double *loci;
  unsigned int total_markers;
  read_chrom(filenames[0],infile_chrom,&chrom_index,&chrom_no_loci,&chrom_length,&chrom_recom_rate,&no_chrom,&loci,&total_markers);

#ifdef DEBUG
  for(int i = 0; i < total_markers; i++)
    {
      printf("%lf\t",loci[i]);
    }
#endif
  chrom_data *head, *current;
  head = data_linked_list_creation(chrom_index, chrom_no_loci, chrom_length, chrom_recom_rate, no_chrom, loci);
  current = head;
  unsigned int n_chromosomes = 0;
  while(current != NULL)
    {
#ifdef DEBUG
      printf("chromosome number:%u\tchromosome length:%lf\trate_of_recom:%lf\tnumber_of_loci:%u\n",current->chrom_id,current->chrom_length, current->chrom_recom_rate, current->n_loci);
      for(int i = 0; i < current->n_loci; i++)
      	{
      	  printf("Marker[%d]:%lf\n",i+1,current->markers[i]);
      	}
#endif
      n_chromosomes++;
      current = current->next;
    }

  /* printf("\nNo. of chromosomes:%u\tTotal no.of markers:%u\n",n_chromosomes,total_markers); */
  char population_name[] = {'A','B','\0'};
  read_pop(filenames[1],infile_A,head,population_name[0]);

#ifdef DEBUG
  current = head;
  while(current != NULL)
    {
      g_hash_table_foreach(current->hashA, (GHFunc)iterator, "Haplotype Sequence: %u has Frequency %lf \n");
      printf("\n\n");
      current = current->next;
    }
#endif

  read_pop(filenames[2],infile_B,head,population_name[1]);

#ifdef DEBUG
  current = head;
  while(current != NULL)
    {
      g_hash_table_foreach(current->hashB, (GHFunc)iterator, "Haplotype Sequence: %u has Frequency %lf \n");
      printf("\n\n");
      current = current->next;
    }
#endif
  
  unsigned int individual_count;
  /* char *names_individuals; */
  /* char individual_name[] */
  read_multiple_indv(filenames[3], infile_indv, head, &individual_count/* , &names_individuals */);

  /* printf("%s\n",names_individuals); */
  
#ifdef DEBUG
  current = head;
  while(current != NULL)
    {
      printf("********* chromosome ********\n\n");
      for(int j = 0; j < individual_count; j++)
	{
	  printf("%u\t%u\n",current->haplotype_1[j], current->haplotype_2[j]);
	}
      current = current->next;
    }
#endif


  double *indv_model_prob[individual_count];
  for(int i = 0; i < individual_count; i++)
    {
      indv_model_prob[i] = (double*) calloc(NO_MODELS, sizeof(double));
      assert(indv_model_prob[i] != NULL);
    }

  double *final_prob[individual_count];
  for(int i = 0; i < individual_count; i++)
    {
      final_prob[i] = (double*) calloc(NO_MODELS, sizeof(double));
      assert(final_prob[i] != NULL);
    }
  

  current = head;
  while(current != NULL)
    {
      update_model_indv_prob(current,indv_model_prob,individual_count);
      for(int i = 0; i < individual_count; i++)
	{
	  for(int j = 0; j < NO_MODELS; j++)
	    {
	      final_prob[i][j] = final_prob[i][j] + indv_model_prob[i][j];
#ifdef DEBUG
	      printf("%lf\t",indv_model_prob[i][j]);
#endif
	    }
	  /* printf("\n"); */
	}
      current = current->next;
      /* printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); */
    }


  /* char row_probabilities[individual_count][BUFFER_SIZE]; */
  /* char new[individual_count][BUFFER_SIZE]; */

  char *row_probabilities[individual_count];
  char *new[individual_count];

  for(int i = 0; i < individual_count; i++)
    {
      row_probabilities[i] = malloc(BUFFER_SIZE);
      new[i] = malloc(BUFFER_SIZE);
      assert(row_probabilities[i] != NULL);
      assert(new[i] != NULL);
    }
  


  /* double log_posterior_denominator[individual_count], scale[individual_count]; */
  double *log_posterior_denominator, *scale;
  log_posterior_denominator = malloc((individual_count) * sizeof(double));
  scale = malloc((individual_count) * sizeof(double));
  assert(log_posterior_denominator != NULL);
  assert(scale != NULL);
  
  for(int i = 0; i < individual_count; i++)
    {
      scale[i] = final_prob[i][0];
      for(int j = 1; j < NO_MODELS; j++)
      	{
	  scale[i] = fmax(final_prob[i][j],scale[i]);
	}
      scale[i] = abs(scale[i])-0.5;
    }



    for(int i = 0; i < individual_count; i++)
    {

      for(int j = 0; j < NO_MODELS; j++)
	{
      log_posterior_denominator[i] = log(exp(final_prob[i][0] + scale[i]) + exp(final_prob[i][1] + scale[i]) + exp(final_prob[i][2] + scale[i]) + exp(final_prob[i][3] + scale[i]) + exp(final_prob[i][4] + scale[i]) + exp(final_prob[i][5] + scale[i])) - scale[i];


	}

    }

    /* printf("The log posterior probabilities ~~~~~~~~~~~~~~~~~~~\n\n"); */
  
  double *log_posterior[individual_count];
  double *posterior_prob[individual_count];
  for(int i = 0; i < individual_count; i++)
    {
      log_posterior[i] = (double*) calloc(NO_MODELS, sizeof(double));
      posterior_prob[i] = (double*) calloc(NO_MODELS, sizeof(double));
      assert(log_posterior[i] != NULL);
      assert(posterior_prob[i] != NULL);
    }

  
  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)
	{
	  log_posterior[i][j] = final_prob[i][j] - log_posterior_denominator[i];
	  posterior_prob[i][j] = exp(log_posterior[i][j]);
	  /* printf("%15lf\t",log_posterior[i][j]); */
	}
      /* printf("\n"); */
    }


  GHashTable *hash_model = g_hash_table_new(g_int_hash, g_int_equal);
  char *model_name_array;
  model_name_array = malloc((NO_MODELS+1) * sizeof(char));
  assert(model_name_array != NULL);
  strcpy(model_name_array,"abcdef");
  int *key;
  key = malloc(NO_MODELS * sizeof(int));
  assert(key != NULL);
  for(int i = 0; i < NO_MODELS; i++)
    {                                   
      key[i] = i;
      g_hash_table_insert(hash_model,&key[i],&model_name_array[i]);
    }                                   
  /* char best_model_post_prob[individual_count+1], second_best_model_post_prob[individual_count+1]; */


  char *best_model_post_prob = NULL;
  best_model_post_prob = malloc(individual_count+1);
  assert(best_model_post_prob != NULL);
 
  char *second_best_model_post_prob = NULL;
  second_best_model_post_prob = malloc(individual_count+1);
  assert(second_best_model_post_prob != NULL);
  
  double *BayesFactor;
  BayesFactor = malloc((individual_count) * sizeof(double));
  assert(BayesFactor != NULL);
  
  for(int i = 0; i < individual_count; i++)
    {
#ifdef DEBUG
      for(int j = 0; j < NO_MODELS; j++)
	{
	  printf("%15lf\t",final_prob[i][j]);
	}
#endif
      largest(posterior_prob[i],NO_MODELS,hash_model,&best_model_post_prob[i],
	      &second_best_model_post_prob[i],&BayesFactor[i]);
      /* printf("\n"); */
    }



  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)

	{
	  /* #ifdef DEBUG */
	  if(j == 0){
	    sprintf(row_probabilities[i],"%15lf\t",final_prob[i][j]);
	  }
	  else{
	    sprintf(new[i],"%15lf\t",final_prob[i][j]);
	    strcat(row_probabilities[i],new[i]);
	  }
	}
      /* strcat(row_probabilities[i],"\n"); */
      /* #endif */
    }


  
  
  /* char append_posterior[individual_count][BUFFER_SIZE]; */
  /* char model_names_BF[individual_count][BUFFER_SIZE]; */

  char *append_posterior[individual_count];
  char *model_names_BF[individual_count];

  for(int i = 0; i < individual_count; i++)
    {
      append_posterior[i] = malloc(BUFFER_SIZE);
      model_names_BF[i] = malloc(BUFFER_SIZE);
      assert(append_posterior[i] != NULL);
      assert(model_names_BF[i] != NULL);
    }

  for(int i = 0; i < individual_count; i++)
    {
      for(int j = 0; j < NO_MODELS; j++)
  	{
	  sprintf(append_posterior[i],"%15lf\t",posterior_prob[i][j]);
	  strcat(row_probabilities[i],append_posterior[i]);
  	}
      sprintf(model_names_BF[i],"\t%8c\t%8c\t%15lf",best_model_post_prob[i],
	      second_best_model_post_prob[i], BayesFactor[i]);
      strcat(row_probabilities[i],model_names_BF[i]);
      strcat(row_probabilities[i],"\n");
    }


  char header[BUFFER_SIZE] = "log_prob_A";
  strcat(header,"\tlog_prob_B");
  strcat(header,"\tlog_prob_C");
  strcat(header,"\tlog_prob_D");
  strcat(header,"\tlog_prob_E");
  strcat(header,"\tlog_prob_F");

  strcat(header,"\tpost_prob_A");
  strcat(header,"\tpost_prob_B");
  strcat(header,"\tpost_prob_C");
  strcat(header,"\tpost_prob_D");
  strcat(header,"\tpost_prob_E");
  strcat(header,"\tpost_prob_F");
  strcat(header,"\tbest_model");
  strcat(header,"\tnext_best");
  strcat(header,"\tBayes_Factor");
  strcat(header,"\n");
  /* printf("%s\n", header); */
  
  gsize bytes_written;
  outfile = g_io_channel_new_file(argv[9],"w",&error);
  g_io_channel_write_chars(outfile,header,strlen(header),&bytes_written,&error);
  for(int i = 0; i < individual_count; i++)
    {

      g_io_channel_write_chars(outfile,row_probabilities[i],strlen(row_probabilities[i]),&bytes_written,&error);
    }
  
  for(int i = 0; i < individual_count; i++)
    {
      free(indv_model_prob[i]);
      free(final_prob[i]);
      free(log_posterior[i]);
      free(posterior_prob[i]);
    }

  free(log_posterior_denominator);
  free(scale);


  free(model_name_array);
  free(key);                                                                                    
  g_hash_table_destroy(hash_model);


  for(int i = 0; i < individual_count; i++)
    {
      free(row_probabilities[i]);
      free(new[i]);
      free(append_posterior[i]);
      free(model_names_BF[i]);
    }

  free(BayesFactor);
  free(best_model_post_prob);
  free(second_best_model_post_prob);

  
  delete_chrom_data(head);
  free(chrom_index);
  free(chrom_no_loci);
  free(chrom_length);
  free(chrom_recom_rate);
  free(loci);


  g_io_channel_unref (infile_chrom);
  g_io_channel_unref (infile_A);
  g_io_channel_unref (infile_B);
  g_io_channel_unref (infile_indv);
  g_io_channel_unref (outfile);
 
  for(int i = 0; i < 4; i++)
    {
      free(filenames[i]);
    }

  
  return 0;
  
}
