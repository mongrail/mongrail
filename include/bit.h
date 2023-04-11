#define MAXLOCI 10
#define MAXHAPS 1024

struct indiv
{
  unsigned int genotype1;
  unsigned int genotype2;
  unsigned int* compHaps;
  int numHaps;
};


