#include <TChain.h>
#include <TChainElement.h>

void printListOfTChainElements(TChain *chain){
  TObjArray *fileElements=chain->GetListOfFiles();
  int nFiles = fileElements->GetEntries();                                                                          
  TIter next(fileElements);
  TChainElement *chEl=0;
  for( int entry=0; entry < nFiles; entry++ ) {
    chEl=(TChainElement*)next();
    printf("%s\n",chEl->GetTitle());
  }
  printf("DEBUG\t\t: %d files in the chain\n",nFiles);
}

void red_Data_BDT(){

  const std::string outfile   = "PhiMuMu_Data_ForFit.root";
  const std::string treename  = "tree";

  TCut b       = "Bmass>4.9 && Bmass<5.8";
  TCut t       = "(JpsiTriggers==1 || PsiPTriggers==1 || LMNTTriggers==1) && (PsiPTriggersdr0<100 || PsiPTriggersdr1<100 || JpsiTriggersdr0<100 || JpsiTriggersdr1<100 || LMNTTriggersdr0<100 || LMNTTriggersdr1<100)";
  TCut t1      = "mtrkqual==1";
  TCut t2      = "ptrkqual==1";
  TCut bl      = "Blxysig>0 && Blxysig<9999";
  

  //resonance rejection
  TCut r1      = "(Mumumass<3.097-5.5*Mumumasserr) || (Mumumass>3.097+3.5*Mumumasserr)";
  TCut r2      = "(Mumumass<3.686-3.5*Mumumasserr) || (Mumumass>3.686+3.5*Mumumasserr)";

  TCut phi     = "Phimass>1.01 && Phimass<1.03";

  TCut trk1 = "Kmtrkdcasigbs>0.8";
  TCut trk2 = "Kptrkdcasigbs>0.8";

  TCut trk  = trk1 && trk2;

  //anti-radiation 
  TCut ar1 = "fabs(Bmass-Mumumass-2.270)>0.120";
  TCut ar2 = "fabs(Bmass-Mumumass-1.681)>0.090";
  
  TCut bdt     = "Bdt>0.14"; 

  TCut mycuts = b && t && t1 && t2 && bl && r1 && r2 && phi && ar1 && ar2 && bdt && trk;

  TChain *ch_sig = new TChain("tree");
  ch_sig->Add("./data/sel_Combine_2016_Mini_Presel_data_cut_bdt-s01.root");

  TTree *inTree = ch_sig;

  if( !inTree ) std::cout << "tree " << treename << " does not exist" << std::endl;
  std::cout << " @@@  entries in sig tree: " << inTree->GetEntries() << std::endl;

inTree->SetBranchStatus("*", 1);

double n_pre = inTree->GetEntries();
double n_pre_wcut = inTree->GetEntries(mycuts);

std::cout << "# of events in original tuple = " << n_pre << endl;
std::cout << "# of events in original tuple after cuts = " << n_pre_wcut << endl;

std::cout << "Copying tree: " << treename.c_str()<< endl;

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
cout << "File created with name-> PhiMuMu_Data_ForFit.root"<< endl;
cout<< "Using the file-> sel_Combine_2016_Mini_Presel_data_cut_bdt-s01.root"<< endl; 


}
