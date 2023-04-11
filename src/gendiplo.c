#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include"bit.h"

/* gcc gendiplo.c -o gendiplo -lm */

void prn_binary(unsigned int x)
{
for (int i = 0; i < 32 && i < MAXLOCI; ++i) {
  if (x >> i & 0x1) putchar('1');
  else putchar('0');
}
}

/* get number of polymorphic sites and number of haplotypes */
unsigned int count_haplotypes(unsigned int hap1, unsigned int hap2)
{
  unsigned int hap1hap2;
  unsigned int no1bits=0;
  hap1hap2 = hap1^hap2; /* 1 if 1^0 or 0^1, 0 otherwise */
  for (int i = 0; i < 32 && i < MAXLOCI; ++i)
    if (hap1hap2 >> i & 0x1) no1bits++; /* if ith bit of hap1hap2 = 1 is polymorphic */
  return pow(2,no1bits); /* no haplotypes = 2^no1bits */
}

/* find all unique haplotypes compatible with genotypes */
void compatible_haps(unsigned int *hapvec, unsigned int genotype1, unsigned int genotype2)
{
  long int g1g2 = genotype1 ^ genotype2;
  int hapNo=1;
  unsigned int mask=1;
  for(int j=0; j<MAXLOCI; j++)
    {
      if((mask << j) & g1g2)
	{
	  hapNo = hapNo*2;
	  for(int k=(hapNo/2);k<hapNo;k++)
	    hapvec[k]=hapvec[k-(hapNo/2)];
	  for(int i=0;i<(hapNo/2);i++)
	    {
	      int allele=0; int pos=0;
	      while(pos<hapNo && allele < 2)
		{
		  if(hapvec[i] == hapvec[pos])
		    {
		      if(allele == 1)
			  hapvec[pos] = (mask << j) ^ hapvec[pos];
		      allele++;
		    }
		  pos++;
		}
	    }
	}
      else
	{
	  if((mask << j) & genotype1)
	    for(int i=0;i<hapNo;i++)
	      hapvec[i] = (mask << j) ^ hapvec[i];
	}
    }
}

/* arrange haplotypes into successive pairs (hap1,hap2) of diplotypes: even = hap1, odd = hap2 */
void sortDiplotypes(struct indiv sampleInd)
{
  int totalDiplos=0;
  int j;
  int i=0;
  if(sampleInd.numHaps > 2)
    {
      while(i < sampleInd.numHaps)
	{
	  j=1;
	  while(j < sampleInd.numHaps)
	    {
	      if((sampleInd.compHaps[i] ^ sampleInd.compHaps[j]) ==
		 (sampleInd.genotype1 ^ sampleInd.genotype2))
		{
		  totalDiplos+=1;
		  unsigned int temp = sampleInd.compHaps[i+1];
		  sampleInd.compHaps[i+1] = sampleInd.compHaps[j];
		  sampleInd.compHaps[j] = temp;
		  i = i+2;
		  break;
		}
	      else
		j++;
	    }
	}
    }
}


int main(int argc, char* argv[])
{
  struct indiv sampleInd;
  int bin_or_dec;
  if(argc == 4)
    {
      bin_or_dec = atoi(argv[1]);
      sampleInd.genotype1 = atoi(argv[2]);
      sampleInd.genotype2 = atoi(argv[3]);
      sampleInd.compHaps = calloc(MAXHAPS,sizeof(unsigned int));
      sampleInd.numHaps = count_haplotypes(sampleInd.genotype1, sampleInd.genotype2);
      compatible_haps(sampleInd.compHaps, sampleInd.genotype1, sampleInd.genotype2);
      sortDiplotypes(sampleInd);
      int numDiplos = sampleInd.numHaps/2;
      
      for(int i=0; i<sampleInd.numHaps; i = i+2)
	if(bin_or_dec == 1)
	  printf("%d %d\n",sampleInd.compHaps[i],sampleInd.compHaps[i+1]);
	else
	  {
	    prn_binary(sampleInd.compHaps[i]);
	    printf(" ");
	    prn_binary(sampleInd.compHaps[i+1]);
	    printf("\n");
	  }
    }
  else
    printf("Usage: gendiplo <bin_or_dec hap1 hap2>\n hap1,hap2 are unsigned integers. bin_or_dec = 0,1\n");
}

