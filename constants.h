//**************//
// Important   //
//  Constants //
//***********//

// Number of iterations to wait for queue to finish
#define MAX_TIME_WAIT 10000
#define MAX_BRANCH_WAIT 10000

// Safe Mode Options
#define SKIPMOPAC 0
#define SKIPDFT 1
#define SKIPGSM 0
#define DO_NOT_WAIT 0
//SKIPSUBMIT is for multistep runs
#define SKIPSUBMIT 1
#define SKIPFIRST 0
#define GEN_SPECIFIC 0

// Screening Thresholds 
#define DFT_EMAX 20     //vs. starting int
#define GSM_EMAX_TS 20  // Ea max = DFT_EMAX + GSM_EMAX
#define MINTSE 10
#define EVSORIGTOL 4.5

//MOPAC/DFT Settings
#define NOMOOPT 1
// use 1 for ASE, 2 for g09
#define USE_ASE_GAUSSIAN 0

//note: use_ase_gaussian ignores all the following settings
#define UNRESTRICTED 1
#define DFT_EXCHANGE "B3LYP"
#define DFT_CORRELATION "NONE"
#define B631GSS 0
#define B631GS  1
#define LANL2DZ 0
#define BASISMIX 0
#define DFT_OPT_CYCLES 285
#define DFT_CHARGE 0
#define DFT_SPIN 1

