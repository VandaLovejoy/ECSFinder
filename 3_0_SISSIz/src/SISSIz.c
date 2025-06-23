/*

  SISSIz -- Dinucleotide controlled random alignments and RNA gene prediction

  Copyright (c) Tanja Gesell & Stefan Washietl, 2008

  See file COPYING for licence.

*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>

#include "dkiss.h"
#include "global.h"
#include "treefile.h"
#include "models.h"
#include "evolve.h"
#include "mutate.h"


#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
#include "alifold.h"
#include "constraints_SHAPE.h"

#include "rnaz_utils.h"
#include "expfit.h"
#include "bionj_custom.h"
#include "cmdline_sissiz.h"
#include "treeML.h"

/* Function prototypes */

TTree* treeFromString(char* treeStrin);
void tree2aln(TTree* tree, struct aln *alignment[]);
void countFreqsDi(const struct aln *alignment[], double freqs[][4]);
void countFreqsMono(const struct aln *alignment[], double freqs[]);
double* localPairID(const struct aln *alignment[]);
double* getCats(double* ids, double level,int L);
double* getCatsFreqs(double* cats, double* ids, double level, int L);

void printAlnClustal(FILE *out, const struct aln* AS[], int printU);
void printAlnMAF(FILE *out, const struct aln* AS[],int printU);
void reintroduceGaps(const struct aln* origAln[], struct aln* sampledAln[]);
double avg(double* data, int N);
double stddev(double* data, int N);
void usage(void);
void help(void);
void version(void);
int compare (const void * a, const void * b);

/* Global variables */
int ntMap[256];

/* Global variables of treefile.c */
extern TNode *avail;
extern long usedAvail;
extern long usedMalloc;



int main(int argc, char *argv[]){

  int i,s,j,k,counter,currRun,runNum,L,N,countTooDiverged;
  double p,d,identity,sum,sumSquared,parA,parB,maxP,difference,adjustment,MFE,mean,stdv;
  double *sampledMFEs;
  double *sampledIDs;
  float sci, gc_content;
    
  double *ids, *targetIDs, *tmpIDs, *pData, *dData, *cats, *catsFreqs, *currCatsFreqs;
  double *Qij=NULL;
  double **QT2ij=NULL;
  double freqsMono[4],freqsDi[4][4];
  double zScore;
 

  int tstvModel=0;
  double kappaPar=4.0;

  double sumFreqsMono[4],sumFreqsDi[4][4];
  double currFreqsMono[4],currFreqsDi[4][4];
  
  double *distMatrix; /* matrix flat in one row */

  char c1,c2;
  char treeString[10000];
  char* treeStringML;
  char *seq1, *seq2, *tmpSeq, *structure;
  
  char** names; 
  
  int seed1=0, seed2=0; //Seed variables for random number generator 
  char** shape_files;
  int with_shapes=0, tmp_number, *shape_file_association;
  unsigned int longest_string;
  char* tmp_string;
  vrna_fold_compound_t *vc;
  char* shapeMethod;
  vrna_md_t md;
  set_model_details(&md);
  
  
  char inputFileName[1024]="STDIN";
  char modelMode[256]="di";
  
  TTree *tree;
  FILE *inputFile=stdin;
  FILE *outputFile=stdout;
  FILE *rateFile;
  FILE *treeFile; /*TMG treeoption 2008*/

	
  int (*readFunction)(FILE *clust,struct aln *alignedSeqs[]);
  struct gengetopt_args_info args;
  struct aln *pair[MAX_NUM_NAMES];     
  struct aln *sampledAln[MAX_NUM_NAMES];     
  struct aln *inputAln[MAX_NUM_NAMES];     
  const char *tmpAln[MAX_NUM_NAMES]; /* simple format for alifold */

  /* Command line arguments and other user defined settings */

  int mockSize=-1;
  int wSSR=1;
  int numSamples=100;
  int simulateOnly=0;
  int verbose=0;
  int mono=0;
  int printRates=0;
  int printTree=0;	/*TMG treeoption 2008*/
  int printU=0;
  int outputMAF=0;
  int regressionSampleSize=10;
  double regressionPointNum=26;
  double dValues[26]={0.0001,0.001,0.01,0.02,0.05,0.08,
                      0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,
                      0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,2,2.5,3};
  double diffMaxP=0.01;

  double precision=1.0;

  double squaredDiff;

  int catGamma=4;

  clock_t beginModelTime,beginSimulationTime;	
  double modelTime,simulationTime;

  beginModelTime=clock();

  for (i=0;i<256;i++) ntMap[i]=0;

  ntMap['A']=ntMap['a']=0;
  ntMap['C']=ntMap['c']=1;
  ntMap['G']=ntMap['g']=2;
  ntMap['T']=ntMap['t']=3;
  ntMap['U']=ntMap['u']=3;
 
  /* Read command line arguments */
  
  if (cmdline_parser (argc, argv, &args) != 0){
    usage();
    exit(EXIT_FAILURE);
  }

  if (args.inputs_num>=1){
    inputFile = fopen(args.inputs[0], "r"); 
    if (inputFile == NULL){
      fprintf(stderr, "ERROR: Can't open input file %s\n", args.inputs[0]);
      exit(1);
    }
    strcpy(inputFileName,args.inputs[0]);
  }
  
  if (args.flanks_given){
    mockSize=args.flanks_arg;
  }
  
  if (args.gamma_given){
    catGamma=args.gamma_arg;
  }

  if (args.rna_given){
    printU=1;
  }

  if (args.dna_given){
    printU=0;
  }

  if (args.maf_given){
    outputMAF=1;
  }
    


  if (args.num_samples_given){
    numSamples=args.num_samples_arg;
  }

  if (args.num_samples_regression_given){
    regressionSampleSize=args.num_samples_regression_arg;
  }

  if (args.precision_given){
    precision=args.precision_arg;
  }

  if (args.outfile_given){
    outputFile = fopen(args.outfile_arg, "w");
    if (outputFile == NULL){
      fprintf(stderr, "ERROR: Can't open output file %s\n", args.outfile_arg);
      exit(1);
    }
  }

  if (args.nossr_given){
    wSSR=0;
  }

  if (args.mono_given){
    mono=1;
    strcpy(modelMode,"mono");
  }
	
   /* TMG alifold 2011*/
  if (args.oldAliEn_given){
		oldAliEn=1;
	}
	
   /* TMG alifold 2011*/	
   if (args.ribo_given){
	    RibosumFile = NULL;
  		 md.ribo     = ribo = 1;
	}
	
  
  if (args.verbose_given){
    verbose=1;
  }

  if (args.simulate_given){
    simulateOnly=1;
    if (!args.num_samples_given){
      numSamples=1;
    }
  }

  if (args.help_given){
    help();
    exit(EXIT_SUCCESS);
  }

  if (args.version_given){
    version();
    exit(EXIT_SUCCESS);
  }

  if (args.print_rates_given){
    printRates=1;
    rateFile=fopen("rates.dat","w");
    if (rateFile==NULL){
      fprintf(stderr,"Could not write file with site specific rates");
      exit(1);
    }
  }

  /*TMG treeoption 2008*/	
  if (args.print_tree_given){    
		printTree=1;
	    treeFile=fopen("aln.tree","w");
		if (treeFile==NULL){
			fprintf(stderr,"Could not write tree output.");
			exit(1);
		} 
	}
	
	
 if (args.tstv_given){
    tstvModel=1;
    if (args.kappa_given){
      kappaPar=args.kappa_arg;
    } else {
      kappaPar=-99.0; /* Estimate from data */
    }
  } else {
    if (args.kappa_given){
      tstvModel=1;
      kappaPar=args.kappa_arg;
    } else {
      kappaPar=1.0;
    }
  }



  /* Read input file */

  switch(checkFormat(inputFile)){
  case CLUSTAL:
    readFunction=&read_clustalz;
    break;
  case MAF:
    readFunction=&read_maf;
    break;
  case 0:
    nrerror("ERROR: Unknown alignment file format. Use Clustal W or MAF format.\n");
  }

  readFunction(inputFile, inputAln);

  for (i=0; inputAln[i]!=NULL; i++){
    tmpSeq=inputAln[i]->seq;
    j=0;
    while (tmpSeq[j]){
      tmpSeq[j]=toupper(tmpSeq[j]);
      j++;
    }
  }

  /* printAlnClustal(outputFile,(const struct aln**)inputAln); */

  /* Calculate important alignment characteristics*/
      
  L=strlen(inputAln[0]->seq); /* Length of alignment without flanks*/

  for (N=0; inputAln[N]!=NULL; N++); /* Number of sequences */

  
  if (N<=1){
    fprintf(stderr,"There must be at least two sequences in the alignment.\n");
    exit(1);
  }

  countFreqsMono((const struct aln**)inputAln, (double *) freqsMono); /* Mononucleotide frequencies */
  countFreqsDi((const struct aln**)inputAln, freqsDi); /* Dinucleotide frequencies */
  
  ids=localPairID((const struct aln**)inputAln); /* Local pairwise identities */

  identity=avg(ids,L); /* Mean pairwise identity of whole alignment */

  targetIDs=(double*)space(sizeof(double)*L); /* Remember original local pairwise ids */



  /* Add "mock sites" */

  if (mockSize<0){
    mockSize=L*1.5; /* Default if no size is given at the command line */
  }

  for (i=0;i<L;i++) targetIDs[i]=ids[i];

  /*cats=getCats(targetIDs,0.01,L);

    catsFreqs=getCatsFreqs(cats, ids, 0.01, L);

    i=0;
    while (catsFreqs[i]>=0.0){
    printf("%.2f ",catsFreqs[i++]);
    }
    printf("\n");
    
  */
  
  free(ids);
  ids=(double*)space(sizeof(double)*(L+mockSize));

  for (i=0;i<L;i++) ids[i]=targetIDs[i]; 
  for (i=L;i<L+mockSize;i++) ids[i]=identity; /* Use average ids for mock-sites*/

  for (i=0;i<L+mockSize;i++) ids[i]=1-ids[i]; /* Turn ids into distances */

  /*if (!simulateOnly || (simulateOnly && verbose)){*/

  if (verbose){
  
    fprintf(outputFile,"# Input file: %s\n",inputFileName);
    fprintf(outputFile,"# Number of sequences: %u\n",N);
    fprintf(outputFile,"# Length: %u\n", L);
    fprintf(outputFile,"# Nucleotide model: %s\n", modelMode);
    /*fprintf(outputFile,"# With site specific sites: %u\n", wSSR);*/
    fprintf(outputFile,"# Transition/Transversion model: %s\n", (tstvModel==1)?("Yes"):("No"));
	/*TMG alifold 2011*/   
	 if (oldAliEn) { 
		 fprintf(outputFile, "# Use old alifold energies (with gaps)\n"); 
	 } else {
		 fprintf(outputFile, "# Use new alifold energies\n"); 
	 }	 
	 if (ribo) fprintf(outputFile, "# Use Ribosum matrices \n");  
	      
	  if (tstvModel){
      if (kappaPar<=0){
        fprintf(outputFile,"# Kappa parameter: Estimate\n");
      } else {
        fprintf(outputFile,"# Kappa parameter: %.2f\n",kappaPar);
      }
    }
    
    fprintf(outputFile,"# Flanking sites: %u\n",mockSize);
    /*fprintf(outputFile,"# No. of samples to simulate: %u\n",numSamples);*/
    fprintf(outputFile,"# No. of samples for regression: %u\n\n",regressionSampleSize);
    fprintf(outputFile,"# Mean pairwise identity: %.4f\n\n",identity);
    fprintf(outputFile,"# Mononucleotide content (ACGT): ");

    for (i=0;i<4;i++) fprintf(outputFile,"%.4f ",freqsMono[i]);
    fprintf(outputFile,"\n\n");
    fprintf(outputFile,"# Dinucleotide content (ACGT x ACGT):\n\n");
    for (i=0;i<4;i++){
      fprintf(outputFile,"# ");
      for (j=0;j<4;j++){
        fprintf(outputFile,"%.4f ",freqsDi[i][j]);
      }
      fprintf(outputFile,"\n");
    }
    fprintf(outputFile,"\n");
  }


  if (tstvModel && kappaPar<0){

    if (verbose){
      fprintf(outputFile,"# Estimating transition/transversion parameter: ");
    }

    treeML((const struct aln**)inputAln,catGamma,&treeStringML,&kappaPar);

    printf("TREE: %s", treeStringML);

    if (verbose){
      fprintf(outputFile,"kappa = %.2f\n",kappaPar);
    }
  }


  /* Set global variables for RNA package */

  do_backtrack=1; 
  dangles=2;

  /* Set global variables for SISSI */
  
  randomgeneratordkiss=1;
  srand (time(NULL));
  

  /* Read seed variables from command line */
  if (args.read_seeds_given) { 
  	seed1=atoi(args.read_seeds_arg[0]);
    if (args.read_seeds_given==2) {
  	    seed2=atoi(args.read_seeds_arg[1]);
        if (seed2==0) {printf("Aborted: Seed value 2 must not be 0"); abort(); }
        start_kiss_pID(seed1,seed2);
    }
    else
  	    {start_kiss(seed1);}
    if (seed1==0) {printf("Aborted:Seed value 1 must not be 0"); abort(); }
  	}
  
  else  { 
  	seed1=(long)rand();
  	seed2=getpid();
 	if (args.seed_with_pID_given)
    		start_kiss_pID(seed1,seed2);
 	 else 
   	 	 start_kiss(seed1);
  }
  

/* Print Seed values */
    if (args.print_seeds_given) {
    FILE *seedFile=fopen("seed.txt","a");
    if (seedFile==NULL){
      fprintf(stderr,"Could not write file with seed values rates");
      exit(1);
    }
  
    if ((args.seed_with_pID_given) || (args.read_seeds_given==2)) 
        fprintf(seedFile,"%i \t% i \n",seed1,seed2);
    else 
        fprintf(seedFile,"%i \n",seed1);
   
    fclose(seedFile);
    }

  nucleotide[0]='A';
  nucleotide[1]='C';
  nucleotide[2]='G';
  nucleotide[3]='U';
  prmat=0;
  prmat2t=0;
 
  if (tstvModel){
    transitiontransverion=1;
    trtv=kappaPar/2;
  } else {
    transitiontransverion=0;
  }

  Rmat[0]=Rmat[1]=Rmat[2]=Rmat[3]=Rmat[4]=Rmat[5]=1.0;
  
  seqlen=L+mockSize;
  numBases=seqlen;

  sitefactor=0;
  verboseWarning=0;
  verbosewaiting=0;
  quiet=1;

  fastqx=1;
  
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      dipi[i*4+j]=freqsDi[i][j];
    }
  }
  
  for (i=0;i<4;i++){
    mpi[i]=freqsMono[i];
  }

  if (!mono){
    di2mono();
    DiMono2Tri();
  }
  
  if (mono){
    Qij=rateRV();
  } else {
    QT2ij=rateT2();
    Pij=computePij(QT2ij);
  }


  /* Sample pairwise alignments for regression */


  dData=(double*)space(sizeof(double)*(regressionPointNum));
  pData=(double*)space(sizeof(double)*(regressionPointNum));

  /* Change A und B for output to be consistent with naming in paper */
  if (verbose) fprintf(outputFile,"\n# Performing regression: p = a * (1 - exp( b * d))\n");

  /* First time without SSRs, second time with SSRs */
  for (runNum=1;runNum<=2;runNum++){
    
    if (verbose && runNum==1) fprintf(outputFile,"# Simulating points (no SSRs)   ");
    if (verbose && runNum==2) fprintf(outputFile,"# Simulating points (with SSRs) ");

    for (counter=0;counter<regressionPointNum;counter++){
    
      d=dValues[counter];

      sprintf(treeString,"(A:%.2f,B:%.2f);",d,d);
      tree=treeFromString(treeString);

      numTaxa=tree->numTips;
      
      sum=0.0;

      for (i=0;i<regressionSampleSize;i++){

        CreateSequences(tree);
        
        if (mono){
          
          EvolveSequencesInd(tree, 1.0, NULL, Qij);

        } else {
        
          EvolveSequencesZ(tree, 1.0, NULL, QT2ij);

        }

        tree2aln(tree,pair);

        tmpIDs=localPairID((const struct aln**)pair);

        p=1-avg(tmpIDs,L+mockSize);
        
        free(tmpIDs);

        freeAln((struct aln**)pair);

        sum+=p;
      }

      
      dData[counter]=d*2;
      pData[counter]=sum/(double)regressionSampleSize;

      if (verbose){
        fprintf(outputFile,".");
        fflush(outputFile);
      }
    
      FreeTree(tree);
      tree=NULL;

      /* Important: also reset these variables from treefile.c */
      avail=NULL;
      usedAvail=0;
      usedMalloc=0;

    }
    
    
    /* Fit the obtained curve and calculate the two parameters of the exponential function:*/
    /* p = B * (1 - exp( A * d)) */

    fitExpFunction(&parA, &parB, counter, pData, dData);

      /* Change A und B for output to be consistent with naming in paper */
    if (verbose) fprintf(outputFile, "  a: %.4f, b: %.4f\n",parB, parA);


    /* In the first run calculate the site specific rates for the second run */
    if (runNum==1){

      /* Starting from local observed differences, SSRs are obtained
         by "correction" using the the formula that has been estimated
         without SSRs before.  Sites than cannot be "corrected" are
         set to the maximum possible divergence level The "missing"
         divergence is then equally distributed over the other
         sites */

      maxP=parB-diffMaxP; 

      countTooDiverged=0;
      difference=0;

      for (i=0;i<L;i++){
        if (ids[i]>maxP){
          difference+=ids[i]-maxP;
          countTooDiverged++;
          ids[i]=maxP;
        }
      }

      adjustment=difference/(double)(L-countTooDiverged);

      for (i=0;i<L+mockSize;i++){
        if (ids[i]<maxP && i<L){
          ids[i]+=adjustment;
        }

        /* Did not see an instance of being again too big after
           adjustment for ungapped case, but in case it happens: */
        
        if (ids[i]>maxP) ids[i]=maxP;

        ids[i]=(1/parA)*log(1-(1/parB)*ids[i]);
                
      }


      /* SSRs are normalized that the sum equals the sequence length */
      sum=0;
      for (i=0;i<L+mockSize;i++){
        if (ids[i]<0.0001) {
          ids[i]=0.0001; /* set zero sites to a very small value */
        }
        sum+=ids[i];
      }

      for (i=0;i<L+mockSize;i++){
        ids[i]=ids[i]/sum*(double)(L+mockSize);
      }


      /* set SISSI global variables and start next round of regression */
      sitefactor=1;
      //scalingfactor=(double *)malloc(seqlen*sizeof(double));

      scalingfactor=ids;
      scalingfactor_sum=(double *)malloc(seqlen*sizeof(double));

      sum=0;
      for(i=0; i<seqlen; i++) {
        scalingfactor_sum[i]=scalingfactor[i]+sum;
        sum=scalingfactor_sum[i];
      }
      
      if (wSSR) continue;
    }
    
    /* If more than two sequences a tree is inferred */
    if (N>2){

      /* Calculate distance matrix */
    
      distMatrix=(double*)space(sizeof(double)*(N*N));
    
      for (i=0;i<N;i++){
        for (j=0;j<N;j++){
          p=0.0;
          if (i==j){
            p=0.0;
          } else {
            seq1=inputAln[i]->seq;
            seq2=inputAln[j]->seq;
            for (k=0;k<L;k++){
              c1=seq1[k];
              c2=seq2[k];
              
              if (c1 == '-' || c2 == '-'){
                p+=targetIDs[k]; /* If gap set to average mean pairwise ID on this site */
              }  else {
                if (c1 == c2){
                  p+=1.0;
                } else {
                  p+=0.0;
                }
              }
            }
            p=1-(p/(double)L);
          }
    
		  /*TMG SW 2008*/
		  /* Set branch length to arbitrary high value if above critical value, to avoid negative logarithm error*/
          if (p>=parB){
          /*    fprintf(stderr,"ERROR: Negative logarithm while generating distance matrix.\n");
            fprintf(stderr,"       Try again with --flanks set to a higher value (currently --flanks=%i)\n",mockSize);
            exit(1);
          }*/
		 /*d=(1/parA)*log(1-(1/parB)*p);*/	  
			  fprintf(stderr, "WARNING: Negative logarithm while gnerating distance matrix.\n");
			  d=0.9;
          } else {
			  d=(1/parA)*log(1-(1/parB)*p);
		  }
		  distMatrix[i*N+j]=d; /* distMatrix is in fact a one dimensonal array */
        }
      }

      if (verbose){
        
        fprintf(outputFile,"\n# Distance matrix\n\n");

        for (i=0;i<N;i++){
          for (j=0;j<N;j++){
            if (j==0) fprintf(outputFile,"# ");
            if (j!=i){
              fprintf(outputFile, "%.4f   ",distMatrix[i*N+j]);
            } else {
              fprintf(outputFile,"  -      ");
            }
          }
          fprintf(outputFile,"\n");
        }
        fprintf(outputFile,"\n");
      }

      treeString[0]='\0';

      /* Perform neighbor joining on the distance matrix */

      treeString[0]='\0';

      bbionj(treeString,(const struct aln**)inputAln,distMatrix);

      if (verbose)   fprintf(outputFile, "\n# BIONJ-Tree\n\n# %s\n\n",treeString);
		
      if (printTree) fprintf(treeFile, "%s\n\n",treeString); /*TMG treeoption 2008*/
	  
	  tree=treeFromString(treeString);

      free(distMatrix);

    } else { /* For pairwise alignments */

      p=(1-identity);

      if (p>=parB){
        fprintf(stderr,"ERROR: Negative logarithm while calculating pairwise distance.\n");
        exit(1);
      }
      
      d=(1/parA)*log(1-(1/parB)*p);

      sprintf(treeString,"(%s:%.f,%s:%.f);",inputAln[0]->name,d/2,inputAln[1]->name,d/2);
    
      tree=treeFromString(treeString);

    }

    if (printRates){
      sum=0.0;
      for (i=0;i<L;i++){
        sum+=ids[i];
      }

      for (i=0;i<L;i++){
        fprintf(rateFile,"%u %.4f\n",i+1,ids[i]/sum*L);
      }
    }

    /* Simulate new alignments along the tree */

    beginSimulationTime=clock();

    numTaxa=tree->numTips;

    /*if (mono){
      Qij=rateRV();
    } else {
      QT2ij=rateT2();
      }*/


    
    if (verbose){
      if (simulateOnly){
        fprintf(outputFile,"# Simulating alignments along the tree...\n");
      } else {
        fprintf(outputFile,"# Simulate and fold alignments... \n\n");
      }
    }

  }

  free(pData);
  free(dData);


  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      sumFreqsDi[i][j]=0.0;
      currFreqsDi[i][j]=0.0;
    }
  }
  
  for (i=0;i<4;i++){
    sumFreqsMono[i]=0.0;
    currFreqsMono[i]=0.0;
  }


  sum=0.0;
  sumSquared=0.0;

  sampledMFEs=space(sizeof(double)*numSamples);
  sampledIDs=space(sizeof(double)*numSamples);

 
  i=0;
    
  //for (i=0;i<numSamples;i++){

  while (i<numSamples){

    CreateSequences(tree);

    if (mono){

      EvolveSequencesInd(tree, 1.0, NULL, Qij);

    } else {

      EvolveSequencesZ(tree, 1.0, NULL, QT2ij);
      
    }

    tree2aln(tree,sampledAln);

    /* Remove mock sites */
    for (j=0;sampledAln[j]!=NULL;j++){
      tmpSeq=sampledAln[j]->seq;
      tmpSeq[L]='\0';
    }

    /* Put gaps back into place */
    reintroduceGaps((const struct aln**)inputAln,(struct aln**)sampledAln);


    countFreqsMono((const struct aln**)sampledAln, currFreqsMono);
    countFreqsDi((const struct aln**)sampledAln, currFreqsDi);

    squaredDiff=0.0;

    for (j=0;j<4;j++) {
      //  printf("%.4f ",currFreqsMono[j]);
      squaredDiff+=(currFreqsMono[j]-freqsMono[j])*(currFreqsMono[j]-freqsMono[j]);
    }

    squaredDiff=sqrt(squaredDiff);

    if (squaredDiff>precision){
      continue;
    } 

    /*for (j=0;j<4;j++)  printf("%.4f ",currFreqsMono[j]);
      printf("%f\n", squaredDiff);*/
    for (j=0;j<4;j++) sumFreqsMono[j]+=currFreqsMono[j];

    for (j=0;j<4;j++){
      for (k=0;k<4;k++){
        sumFreqsDi[j][k]+=currFreqsDi[j][k];
      }
    }

    tmpIDs=localPairID((const struct aln**)sampledAln);

    /*currCatsFreqs=getCatsFreqs(cats, tmpIDs, 0.01, L);

    j=0;

    squaredDiff=0.0;

    j=0;
    while (catsFreqs[j]>=0.0){
      squaredDiff+=(currCatsFreqs[j]-catsFreqs[j])*(currCatsFreqs[j]-catsFreqs[j]);
      j++;
    }

    squaredDiff=sqrt(squaredDiff);

    printf(" %.2f\n",squaredDiff);

    */

    sampledIDs[i]=avg(tmpIDs,L);

    free(tmpIDs);

    if (simulateOnly){

      if (outputMAF){

        for (j=0;sampledAln[j]!=NULL;j++){
          for (k=0;inputAln[k]!=NULL;k++){
            if (strcmp(sampledAln[j]->name,inputAln[k]->name)==0){
              sampledAln[j]->start=inputAln[k]->start;
              sampledAln[j]->length=inputAln[k]->length;
              sampledAln[j]->fullLength=inputAln[k]->fullLength;
              sampledAln[j]->strand=inputAln[k]->strand;
              break;
            }
          }
        }
        printAlnMAF(outputFile,(const struct aln**)sampledAln,printU);
      } else {
        printAlnClustal(outputFile,(const struct aln**)sampledAln,printU);
      }

      
    } else {
      
      structure = (char *) space((unsigned) L+1);
		
      for (j=0;sampledAln[j]!=NULL;j++){
        tmpAln[j]=sampledAln[j]->seq;
      }

      tmpAln[j]=NULL;
	  vc = vrna_fold_compound_comparative((const char **) tmpAln, &md, VRNA_OPTION_DEFAULT);
      MFE = vrna_mfe(vc, structure);

      sum+=MFE;
      sumSquared+=MFE*MFE;

      sampledMFEs[i]=MFE;
      
      if (verbose) fprintf(outputFile, "# %u: %s %.2f\n\n",i, structure, MFE);
      
      free(structure);

    }

    freeAln(sampledAln);

    i++;

  }

  FreeTree(tree);

  for (j=0;j<4;j++) sumFreqsMono[j]/=numSamples;

  for (j=0;j<4;j++){
    for (k=0;k<4;k++){
      sumFreqsDi[j][k]/=numSamples;
    }
  }

  if (verbose){
  
    fprintf(outputFile,"# Sampled mononucleotide content (ACGT): ");
    for (i=0;i<4;i++) fprintf(outputFile,"%.4f ",sumFreqsMono[i]);
    fprintf(outputFile,"\n\n");
    fprintf(outputFile,"# Sampled dinucleotide content (ACGT x ACGT):\n\n");
    for (i=0;i<4;i++){
      fprintf(outputFile,"# ");
      for (j=0;j<4;j++){
        fprintf(outputFile,"%.4f ",sumFreqsDi[i][j]);
      }
      fprintf(outputFile,"\n");
    }
    fprintf(outputFile,"\n");

    /*
      for (i=0;i<4;i++){
      fprintf(outputFile,"# ");
      for (j=0;j<4;j++){
      fprintf(outputFile,"%.4f",freqsMono[i]*freqsMono[j],freqsMono[i]);
      }
      fprintf(outputFile,"\n");
      }
    */
  
    fprintf(outputFile,"# Average mean pairwise identity of samples: %.4f\n",avg(sampledIDs,numSamples));

    if (numSamples>2){
      fprintf(outputFile,"# Std. dev. of mean pairwise identites of samples: %.4f\n\n",stddev(sampledIDs,numSamples));
    }
  }

  if (!simulateOnly){
  	
     /* SHAPE reactivity data */ //Adapted from VIENNA RNA PACKAGE RNAalifold.c (Version 2.3.2 )
    if(args.shape_given) {
    
	    printf("SHAPE FILE GIVEN\n"); 

	    if (verbose)
	      vrna_message_info(stderr, "SHAPE reactivity data correction activated");

	    with_shapes = 1;
	    shape_files = (char **)vrna_alloc(sizeof(char *) * (args.shape_given + 1));
	    shape_file_association = (int *)vrna_alloc(sizeof(int *) * (args.shape_given + 1));

	    /* find longest string in argument list */
	    longest_string = 0;
	    for (s = 0; s < args.shape_given; s++)
	      if (strlen(args.shape_arg[s]) > longest_string)
		longest_string = strlen(args.shape_arg[s]);

	    tmp_string  = (char *)vrna_alloc(sizeof(char) * (longest_string + 1));
	    tmp_number  = 0;

	    for (s = 0; s < args.shape_given; s++) {
	      /* check whether we have int=string style that specifies a SHAPE file for a certain sequence number in the alignment */
	      if (sscanf(args.shape_arg[s], "%d=%s", &tmp_number, tmp_string) == 2) {
		shape_files[s]            = strdup(tmp_string);
		shape_file_association[s] = tmp_number - 1;
	      } else {
		shape_files[s]            = strdup(args.shape_arg[s]);
		shape_file_association[s] = s;
	      }

	      if (verbose) {
		vrna_message_info(stderr,
		                  "Using SHAPE reactivity data provided in file %s for sequence %d",
		                  shape_files[s],
		                  shape_file_association[s] + 1);
	      }
    	  }

    shape_file_association[s] = -1;
    free(tmp_string);
    }
  
    else  if (verbose) {printf("SHAPE FILE NOT GIVEN\n"); }
    
    if (with_shapes) {
      for (s = 0; shape_file_association[s] != -1; s++);

      if (s != N)
        vrna_message_warning("Number of sequences in alignment does not match number of provided SHAPE reactivity data files! ");

      shape_files             = (char **)vrna_realloc(shape_files, (N + 1) * sizeof(char *));
      shape_file_association  = (int *)vrna_realloc(shape_file_association, (N + 1) * sizeof(int));
    }
    if (with_shapes) {
      if (args.shapeMethod_given) {
      	shapeMethod=strdup(args.shapeMethod_arg);
      	}
      else {
      	shapeMethod="D"; 
      	}
      	
      for (i=0;inputAln[i]!=NULL;i++){
      		tmpAln[i]=inputAln[i]->seq;
   	    }
      vc = vrna_fold_compound_comparative((const char **) tmpAln, &md, VRNA_OPTION_DEFAULT);
      
      vrna_constraints_add_SHAPE_ali(vc, \
      shapeMethod, \
      (const char **)shape_files, \
      shape_file_association, \
      verbose, \
      0);  //VRNA_OPTION_MFE | ((pf) ? VRNA_OPTION_PF : 0) //pf ... partition funciton not used in SISSIZ
      
    
    structure = (char *) space((unsigned) L+1);
    MFE = vrna_mfe(vc, structure);
    //MFE = vrna_eval_structure(vc, structure);
   
    //printf("MFE%fMFE",MFE);
    
    }
    
    if (!with_shapes) {
     	for (i=0;inputAln[i]!=NULL;i++){
      		tmpAln[i]=inputAln[i]->seq;
   	}
    	vc = vrna_fold_compound_comparative((const char **) tmpAln, &md, VRNA_OPTION_DEFAULT);
    	structure = (char *) space((unsigned) L+1);
    	MFE = vrna_mfe(vc, structure);
    }
    free(structure);

    mean=avg(sampledMFEs,numSamples);
    stdv=stddev(sampledMFEs,numSamples);
    
    zScore=0.0;

    if (stdv > 0.0){
      zScore=(MFE-mean)/stdv;
    }

    //Calculate SCI 
    if (args.sci_given) {
    float e_mean; 
    sci= MFE; 
    for (i=0;inputAln[i]!=NULL;i++){
        tmpAln[i]=inputAln[i]->seq;
        char *seq_temp = get_ungapped_sequence(tmpAln[i]);
        e_mean += vrna_fold(seq_temp, NULL);
        free(seq_temp);
    }
    e_mean /= i;
    if (e_mean == 0.)
        sci = 0.;
    else
        sci /= e_mean;
    if (sci<0.1) vrna_message_warning("SCI below 0.1"); 
    }
    // Calculate GC content 
    
    float gaps_temp=0; 
    for (i=0;inputAln[i]!=NULL;i++){
        tmpAln[i]=inputAln[i]->seq; 
        for (j=0;j<L;j++) {
            if ( (tmpAln[i][j]=='G') || (tmpAln[i][j]=='C') ) 
                gc_content++; 
            if (tmpAln[i][j]=='-')
                gaps_temp++;
        }
    }
    gc_content/=(N*L-gaps_temp); 

  }
  
  if (!simulateOnly && verbose){
    fprintf(outputFile, "# Consensus MFE native alignment: %.2f\n",MFE);
    fprintf(outputFile, "# Average of sampled alignments: %.2f\n",mean);
    fprintf(outputFile, "# Std. dev. of sampled alignments: %.2f\n",stdv);
    fprintf(outputFile, "# z-score: %.2f\n",zScore);
  }
    
  

  if (!simulateOnly){
    if (args.sci_given) 
        fprintf(outputFile,"sissiz-%s\t%s\t%u\t%u\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\n",modelMode,inputFileName,N,L,identity,avg(sampledIDs,numSamples),stddev(sampledIDs,numSamples),sci,gc_content,MFE,mean,stdv,zScore);
    else
        fprintf(outputFile,"sissiz-%s\t%s\t%u\t%u\t%.4f\t%.4f\t%.4f\tNA\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\n",modelMode,inputFileName,N,L,identity,avg(sampledIDs,numSamples),stddev(sampledIDs,numSamples),gc_content,MFE,mean,stdv,zScore);
  }


  free(ids);
  free(targetIDs);
  free(sampledMFEs);
  //free(sampledIDs);

  if (mono){
    free(Qij);
  } else {
    for (i=0;i<64;i++){
      free(QT2ij[i]);
    }
    free(QT2ij);
  }

  freeAln((struct aln**)inputAln);

  cmdline_parser_free(&args);

  if (verbose){

    modelTime = (double)(clock() - beginModelTime) / CLOCKS_PER_SEC;
    simulationTime=(double)(clock() - beginSimulationTime) / CLOCKS_PER_SEC;

    fprintf(outputFile,"\n# Time used to build model:   %.2fs\n",modelTime-simulationTime);
    fprintf(outputFile,"# Time used to simulate/fold: %.2fs\n",simulationTime);
    fprintf(outputFile,"# Total time used:            %.2fs\n",modelTime);

  }



  exit(EXIT_SUCCESS);

}

TTree* treeFromString(char* treeString){

  FILE * tmpfh;
  TTree* tree;
  int dummy1;
  double dummy2;

  tmpfh = tmpfile();
  fprintf(tmpfh,"%s",treeString);
  rewind(tmpfh);
  tree=NewTree();
  InitTree(tree);
  ReadTree(tmpfh, tree, 0, 0, NULL, &dummy1, &dummy2);
  fclose(tmpfh);
  return tree;
 
}

void tree2aln(TTree* tree, struct aln *alignment[]){

  int i;

  for (i=0; i<tree->numTips; i++) {
    alignment[i]=createAlnEntry(strdup(tree->names[i]), 
                                strdup(tree->tips[i]->sequence),0,0,0,'?');
  
  }
  alignment[i]=NULL;
}


void countFreqsMono(const struct aln *alignment[], double freqs[]){

  int i,k;
  char* currSeq;
  char c;
  unsigned long counter;
  double sum;

  for (i=0;i<4;i++) freqs[i]=0.0;

  counter=0;

  for (i=0; alignment[i]!=NULL; i++){

    currSeq=alignment[i]->seq;

    k=0;
    while ((c=currSeq[k++])!='\0'){
      if (c=='-') continue;
      freqs[ntMap[c]]++;
      counter++;
    }
  }

  for (i=0;i<4;i++) freqs[i]/=(double)counter;

  /*
  for (i=0;i<4;i++){
    printf("%.4f ",freqs[i]);
  }
  printf("\n");
  
  sum=0;

  for (i=0;i<4;i++) sum+=freqs[i];

  printf("Sum: %f\n",sum);
  */

}

void countFreqsDi(const struct aln *alignment[], double freqs[][4]){

  int i,j,k;
  char* currSeq;
  char c1,c2;
  unsigned long counter;
  double sum;
  
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){ 
      freqs[i][j]=0.0;
    }
  }

  counter=0;

  for (i=0; alignment[i]!=NULL; i++){

    currSeq=alignment[i]->seq;

    for (k=0; k<strlen(currSeq)-1;k++){
      c1=currSeq[k];
      c2=currSeq[k+1];

      if (c1!='-' && c2!='-'){
        freqs[ntMap[c1]][ntMap[c2]]+=1;
        counter++;
      }
    }
  }

  for (i=0;i<4;i++){
    for (j=0;j<4;j++){ 
      freqs[i][j]/=(double)counter;
    }
  }


  /*
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      printf("%.4f ",freqs[i][j]); 
    }
    printf("\n");
  }
  printf("\n");
  sum=0;
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){ 
      sum+=freqs[i][j];
    }
  }
  printf("Sum di: %f\n",sum);
  */

}

double* localPairID(const struct aln *alignment[]){

  int i,j,k,length;
  int pairs,matches;

  double* ids;

  length=strlen(alignment[0]->seq);

  ids = (double *) space((sizeof(double))*length);

  for (k=0;k<length;k++){
    
    pairs=0;
    matches=0;
    
    for (i=0;alignment[i]!=NULL;i++){
      for (j=i+1;alignment[j]!=NULL;j++){
        /* Ignore comparisons with gaps */
        if ((alignment[i]->seq[k]!='-') && (alignment[j]->seq[k]!='-')){
          if (alignment[i]->seq[k]==alignment[j]->seq[k]){
            matches++;
          }
          pairs++;
        }
      }
    }

    /* only one non-gap character in column (or gaps only, but this should not happen)*/
    if (pairs==0){
        ids[k]=1;
    } else {
      if (matches>0){
        ids[k]=(double)(matches)/pairs;
      } else {
        ids[k]=0;
      }
    }
  }

  /* for (k=0;k<length;k++){
    printf("%.4f ",ids[k]);
  }

  printf("\n");
  */

  return ids;

}


double avg(double* data, int N){

  int i;
  double sum;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=data[i];
  }
  
  return sum/(double)N;

}

double stddev(double* data, int N){

  int i;
  double sum;
  double mean;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=data[i];
  }

  mean=sum/(double)N;

  sum=0.0;

  for (i=0;i<N;i++){
    sum+=(mean-data[i])*(mean-data[i]);
  }
  
  return sqrt(sum/(double)(N-1));
}


void printAlnClustal(FILE *out, const struct aln* AS[],int printU){

  int i,j,N;
  int nameLength=0;
  int width=60;
  char **paddedNames;
  char* tmpString;
  int start, end;
  int L;

  L=strlen(AS[0]->seq);

  for (i=0;AS[i]!=NULL;i++){
    if (strlen(AS[i]->name)>nameLength){
      nameLength=strlen(AS[i]->name);
    }
  }
  
  for (i=0;AS[i]!=NULL;i++){

    tmpString=AS[i]->seq;

    j=0;
    
    while (tmpString[j]){
      if (!printU){
        if (tmpString[j]=='U'){
          tmpString[j]='T';
        }
      }
      j++;
    }
  }
  
  tmpString=(char*)space(sizeof(char)*(L+1));

  N=i;

  paddedNames=(char**)space(sizeof(char*)*N);

  for (i=0;i<N;i++){
    paddedNames[i]=(char*)space(sizeof(char)*(nameLength+1));
    strcpy(paddedNames[i],AS[i]->name);
    for (j=strlen(AS[i]->name);j<nameLength;j++){
      paddedNames[i][j]=' ';
    }
    paddedNames[i][j]='\0';
  }

  end=0;
  start=0;

  fprintf(out,"CLUSTAL W (SISSIz %s simulation)\n\n\n",PACKAGE_VERSION);

  while (1){
    
    end=start+width;

    if (end>L){
      end=L;
    }

    for (i=0;AS[i]!=NULL;i++){
      fprintf(out, "%s        ",paddedNames[i]);
      strncpy(tmpString,(AS[i]->seq)+start,end-start);
      
      tmpString[end-start]='\0';
      
      fprintf(out, "%s\n",tmpString);

    }
    fprintf(out, "\n\n");

    start=end;

    if (start==L) break;

  }

  for (i=0;AS[i]!=NULL;i++){
    free(paddedNames[i]);
  }

  free(tmpString);
  free(paddedNames);

}


void printAlnMAF(FILE *out, const struct aln* AS[],int printU){

  int i,j,N;
  int L;
  char* tmpString;

  L=strlen(AS[0]->seq);

  fprintf(out, "a score=0\n");

  for (i=0;AS[i]!=NULL;i++){

    tmpString=AS[i]->seq;

    j=0;
    
    while (tmpString[j]){
      if (!printU){
        if (tmpString[j]=='U'){
          tmpString[j]='T';
        }
      }
      j++;
    }
  }
  
  for (i=0;AS[i]!=NULL;i++){
    fprintf(out, "s %s %i %i %c %i %s\n",AS[i]->name,AS[i]->start,AS[i]->length,AS[i]->strand, AS[i]->fullLength, AS[i]->seq);
  }
  fprintf(out, "\n");
  
}



void reintroduceGaps(const struct aln* origAln[], struct aln* sampledAln[]){

  int i,j,k;

  char* tmpSeq;
  char* tmpName;
  char* origSeq;
  char* sampledSeq;

  //printAlnClustal((const struct aln**)origAln);


  for (i=0;origAln[i]!=NULL;i++){
    for (j=0;sampledAln[j]!=NULL;j++){

      //      printf("%s| %s| %u %u\n",origAln[i]->name,sampledAln[j]->name,strlen(origAln[i]->name),strlen(sampledAln[j]->name));


      if (strcmp(origAln[i]->name,sampledAln[j]->name)==0){

        tmpSeq=sampledAln[j]->seq;
        tmpName=sampledAln[j]->name;

        sampledAln[j]->seq=sampledAln[i]->seq;
        sampledAln[j]->name=sampledAln[i]->name;

        sampledAln[i]->seq=tmpSeq;
        sampledAln[i]->name=tmpName;

        origSeq=origAln[i]->seq;
        sampledSeq=sampledAln[i]->seq;

        for (k=0;k<strlen(origSeq);k++){
          if (origSeq[k]=='-'){
            sampledSeq[k]='-';
          }

        }
      }
    }
  }
}


double* getCats(double* ids, double level, int L){

  double currCat, lastCat;
  
  double* sorted;
  double* cats;

  int i,j;

  sorted=(double*)space(sizeof(double)*L);
  cats=(double*)space(sizeof(double)*(L+1));

  for (i=0;i<L;i++) sorted[i]=ids[i];

  qsort(sorted, L, sizeof(double), compare);  
  
  lastCat=sorted[0];
  cats[0]=sorted[0];
  j=1;

  for (i=1;i<L;i++){
    currCat=sorted[i];
    if (currCat-lastCat>level){
      cats[j]=currCat;
      lastCat=currCat;
      j++;
    }
  }

  for (i=0;i<L;i++){
    //printf("%.2f ",sorted[i]);
  }

  cats[j]=-1.0;

  return cats;
    
}


double* getCatsFreqs(double* cats, double* ids, double level, int L){

  int i,j;

  double* freqs;

  freqs=(double*)space(sizeof(double)*(L+1));

  i=0;

  while (cats[i]>=0.0){
    freqs[i]=0.0;
    i++;
  }

  i=0;

  while (cats[i]>=0.0){

    for (j=0;j<L;j++){

      if (fabs(cats[i]-ids[j])<level){
        freqs[i]+=1.0;
      }
    }
    i++;
  }

  freqs[i]=-1;
  
  i=0;

  while (cats[i]>=0.0){
    freqs[i]/=L;
    //  printf("%.2f: %.2f ",cats[i],freqs[i]);
    i++;
  }

  return freqs;

}





int compare (const void * x,const void * y)
{
  if (*(double*)x > *(double*)y)
      return 1;
  else if (*(double*)x < *(double*)y)
      return -1;
   else
      return 0;
}


void usage(void){
  help();
}

void help(void){

  cmdline_parser_print_version ();

  printf("\nUsage: %s [OPTIONS]... [FILE]\n\n", CMDLINE_PARSER_PACKAGE);
  printf("%s\n","  -n, --num-samples    Number of alignments to be sampled");
  printf("%s\n","  -d, --di             Dinucleotide model (default)");
  printf("%s\n","  -i, --mono           Mononucleotide model");
  printf("%s\n","  -s, --simulate       Simulate only (no folding)");
  printf("%s\n","  -v, --verbose        More verbose screen output");
  printf("%s\n","  -t, --tstv           Consider transitions transversion model");
  printf("%s\n","  -k, --kappa          Kappa (default: estimated from data)");
  printf("%s\n","  -p, --precision      Set precision of monunucleotide content");
  printf("%s\n","  --dna, --rna         Print Ts (default) or Us");
  printf("%s\n","  --clustal, --maf     CLUSTAL (default) or MAF output");
  printf("%s\n","  -m, --num-regression Number of sampled points for regression");
  printf("%s\n","  -f, --flanks         Number of flanking 'buffer' sites");
  printf("%s\n","  -b, --print-tree     Print BIONJ-Tree to aln.tree");	  /*TMG treeoption 2008*/
  printf("%s\n","  -a, --oldAliEn       Use old alifold energies");	  /*TMG alifold 2011*/
  printf("%s\n","  -j, --ribo           Ribosome matrices of Alifold are used.");	  /*TMG alifold 2011*/
  printf("%s\n","  -x, --print-rates    Print rates to rates.dat (debugging)");
  printf("%s\n","  --shape              Use SHAPE reactivity data to guide structure predicitons =file1,file2");
  printf("%s\n","  --shapeMethod        Specify the method how to convert SHAPE reactivity data to pseudo energy contributions =D[mX][bY]");
  printf("%s\n","  -y, --seed_with_pID  Seed with process ID. Default: OFF");
  printf("%s\n","  -z, --print_seeds    Print seed values to seed.txt . Default: OFF");
  printf("%s\n","  --read_seeds         Set seed value(s) ={seed1(,seed2)}.");
  printf("%s\n","  --sci                Calculate SCI. Default: OFF");
  printf("%s\n","  -h, --help           Help screen");
  printf("%s\n\n","  -V, --version      Print version");
}

void version(void){
  printf("SISSIz v " PACKAGE_VERSION "\n");
  exit(EXIT_SUCCESS);
}

