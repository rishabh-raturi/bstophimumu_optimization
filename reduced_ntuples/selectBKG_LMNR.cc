#include <TChain.h> 
#include <TChainElement.h>

void printListOfTChainElements(TChain *chain){
  TObjArray *fileElements=chain->GetListOfFiles();
  int nFiles = fileElements->GetEntries();
  //printf("DEBUG\t\t: %d files in the chain\n",nFiles);
  TIter next(fileElements);
  TChainElement *chEl=0;
  for( int entry=0; entry < nFiles; entry++ ) {
    chEl=(TChainElement*)next();
    printf("%s\n",chEl->GetTitle());
  }
  printf("DEBUG\t\t: %d files in the chain\n",nFiles);
}

void selectBKG_LMNR(){ 

  // Nov14: added trigger cuts
  // Major difference:

  //   1. consider 1 < m(dimu) < 2.7 || m(dimu) > 4 GeV // till the time discrepancy in trigger is removed
  //   2. sidebands from data: (3.5, 6.0)*sigma and (3.5, 6.0)*sigma of Bs pdg mass

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // -- define input file name, output file name and cuts to apply
  // -- and also the name of the output file
  // -- anti-radiation cuts not to be applied as they are optimized after BDT cut value is obtained.
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
  const std::string outfile   = "BKG_ForBDT_Scan.root";
  const std::string treename  = "tree";
  
  TCut t = "(JpsiTriggers==1 || PsiPTriggers==1 || LMNTTriggers==1) && (PsiPTriggersdr0<100 || PsiPTriggersdr1<100 || JpsiTriggersdr0<100 || JpsiTriggersdr1<100 || LMNTTriggersdr0<100 || LMNTTriggersdr1<100)";
  TCut di = "(Mumumass>=1.0 && Mumumass<=2.9) || (Mumumass>=4.0 && Mumumass<=4.9)";

  //resonance rejection
  TCut r1      = "(Mumumass<3.096-5.5*Mumumasserr) || (Mumumass>3.096+3.5*Mumumasserr)";
  TCut r2      = "(Mumumass<3.686-3.5*Mumumasserr) || (Mumumass>3.686+3.5*Mumumasserr)";

  //side-band
  TCut sb      = "(Bmass>=5.189724 && Bmass<=5.263589) || (Bmass>=5.470411 && Bmass<=5.544276)";

  TCut phi     = "Phimass>1.01 && Phimass<1.03";

  //track-quality
  TCut t1      = "mtrkqual==1";
  TCut t2      = "ptrkqual==1";

  TCut bl      = "Blxysig>0 && Blxysig<9999";
 

  TCut mycutb  = t  && r1 && r2 && sb && phi && t1 && t2 && bl;

  
// -------------------- Copy the tree -----------------------------------
    

TChain *ch = new TChain("tree");
ch->Add("./data/sel_Combine_2016_Mini_Presel_data_cut_bdt-s01.root");
printListOfTChainElements(ch);

TTree *inTree = ch;

if( !inTree ) std::cout << "tree " << treename << " does not exist" << std::endl;
std::cout << " @@@  entries in bkg tree: " << inTree->GetEntries() << std::endl;


// -- activate the branches you need  
inTree->SetBranchStatus("*", 1);   
      
double n_pre = inTree->GetEntries();
double n_pre_wcut = inTree->GetEntries(mycutb);

std::cout << "# of events in original tuple = " << n_pre << endl;  
std::cout << "# of events in original tuple after cuts = " << n_pre_wcut << endl;  

std::cout << "Copying tree: " << treename.c_str() << endl;
    
TFile* newFile = new TFile(outfile.c_str(),"RECREATE");
///TTree* outTree = inTree->CopyTree( cuts.c_str() );
TTree* outTree = inTree->CopyTree( mycutb );
    
int percentCounter = 1;

for(int i = 0; i < inTree->GetEntries(); ++i){
  
  const int percent = (int)(inTree->GetEntries()/100.0);
    
  if( i == percent*percentCounter ){
    std::cout << percentCounter << " %" << std::endl;
    percentCounter++;
  }
      
 }
    
outTree->Write();
newFile->Save();

double n_post = outTree->GetEntries();
    
std::cout << "# of events in the background tuple = " << n_post << endl;


} 


