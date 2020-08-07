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

void red_JpsiMC_BDT()
{
  const std::string outfile   = "JPsi_MC_ForFit.root";
  const std::string treename  = "tree";

  TCut b       = "Bmass>4.9 && Bmass<5.8";
  
  TCut t       = "(JpsiTriggers==1 ) && ( JpsiTriggersdr0<100 || JpsiTriggersdr1<100 )";
  TCut t1      = "mtrkqual==1";
  TCut t2      = "ptrkqual==1";
  TCut bl      = "Blxysig>0 && Blxysig<9999";

  //Jpsi selection
  TCut jpsi    = "Mumumass<3.097+3.0*Mumumasserr && Mumumass>3.097-5.0*Mumumasserr";


  TCut phi     = "Phimass>1.01 && Phimass<1.03";

  TCut bdt     = "Bdt>0.14"; 

  TCut trk1 = "Kmtrkdcasigbs>0.8";
  TCut trk2 = "Kptrkdcasigbs>0.8";
  
  TCut trk  = trk1 && trk2;

  TCut mycuts = b && t && t1 && t2 && bl && jpsi && phi && bdt && trk;

  TChain *ch_sig = new TChain("tree");
  ch_sig->Add("./sel_BsToJpsiPhi_2016MC_Official_Presel_mc.lite_cut_bdt-s01.root");

  TTree *inTree = ch_sig;

  if( !inTree ) std::cout << "tree " << treename << " does not exist" << std::endl;
  std::cout << " @@@  entries in sig tree: " << inTree->GetEntries() << std::endl;

inTree->SetBranchStatus("*", 1);

double n_pre = inTree->GetEntries();
double n_pre_wcut = inTree->GetEntries(mycuts);

std::cout << "# of events in original tuple = " << n_pre << endl;
std::cout << "# of events in original tuple after cuts = " << n_pre_wcut << endl;

std::cout << "Copying tree: " << treename.c_str() << endl;

TFile* newFile = new TFile(outfile.c_str(),"RECREATE");
TTree* outTree = inTree->CopyTree( mycuts );

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

std::cout << "# of events in the signal tuple = " << n_post << endl;


}

