#include "TTH/StandaloneBDT/interface/BDTClassifier.h"
#include "TTH/CommonClassifier/interface/MEMClassifier.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TGraph.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TTreeFormula.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"

typedef vector<TLorentzVector> LVs;
typedef TLorentzVector LV;

// helper functions to get LorentzVectors of jets
TLorentzVector getLV(float pt,float eta,float phi,float e){
  TLorentzVector v;
  v.SetPtEtaPhiE(pt,eta,phi,e);
  return v;
}

LVs getLVs(uint n, float* pt,float* eta,float* phi,float* e){
  vector<TLorentzVector> vs;
  for(uint i=0; i<n; i++){
    TLorentzVector v= getLV(pt[i],eta[i],phi[i],e[i]);
    vs.push_back(v);
  }
  return vs;
}



class Evaluater{
public:
  Evaluater(const TString infileName, const TString outfileName);
  ~Evaluater();
  void eval();
  
private:
  
  TString infileName;
  TString outfileName;

  
  TTree* inTree;
  TTree* outTree;
TFile* inTreeFile;
 TFile* outTreeFile; 
  
  Float_t BDToutput;
  Float_t MEoutput;
  Float_t p_sig;
  Float_t p_bkg;

  //Integration uncertainties of the probabilities
  Float_t p_err_sig;
  Float_t p_err_bkg;
  //Integration uncertainties of the probabilities
  Float_t n_perm_sig;
  Float_t n_perm_bkg;

 Int_t category; 
 
 Float_t jetPt[100]; 
  Float_t jetPhi[100]; 
 Float_t jetEta[100]; 
 Float_t jetE[100]; 
 Float_t jetCSV[100]; 

  Float_t leptonPt; 
  Float_t leptonPhi; 
  Float_t leptonEta; 
  Float_t leptonE;
  
  Float_t metPt;
  Float_t metPhi;
  
  Float_t weight;
  
  Int_t nJets;
  Int_t evtID;  

  Float_t oldBDT;
  
  BDTClassifier* BDT;
  MEMClassifier* MEM;
  
   std::map<std::string,float> bdtVarMap;

   TStopwatch* timer = new TStopwatch();
   TStopwatch* totaltimer = new TStopwatch();
};

Evaluater::Evaluater(const TString infileName, const TString outfileName){
  
  std::cout<<"files "<<infileName<<" "<<outfileName<<std::endl;
  
  inTreeFile=new TFile(infileName,"READ");
  inTree=(TTree*)inTreeFile->Get("MVATree");

  inTree->SetBranchAddress("Jet_Pt",&jetPt);
  inTree->SetBranchAddress("Jet_Eta",&jetEta);
  inTree->SetBranchAddress("Jet_Phi",&jetPhi);
  inTree->SetBranchAddress("Jet_E",&jetE);
  inTree->SetBranchAddress("Jet_CSV",&jetCSV);

  inTree->SetBranchAddress("N_Jets",&nJets);
  
  inTree->SetBranchAddress("Evt_Pt_PrimaryLepton",&leptonPt);
  inTree->SetBranchAddress("Evt_Eta_PrimaryLepton",&leptonEta);
  inTree->SetBranchAddress("Evt_Phi_PrimaryLepton",&leptonPhi);
  inTree->SetBranchAddress("Evt_E_PrimaryLepton",&leptonE);
  
  inTree->SetBranchAddress("Evt_Pt_MET",&metPt);
  inTree->SetBranchAddress("Evt_Phi_MET",&metPhi);

  inTree->SetBranchAddress("Weight",&weight);
  inTree->SetBranchAddress("Evt_ID",&evtID);
  
  inTree->SetBranchAddress("BDT_v4_output",&oldBDT);
  
  outTreeFile=new TFile(outfileName,"RECREATE");
  outTreeFile->cd();
  outTree=new TTree("MVATree","MVATree");
  outTree->Branch("BDToutput",&BDToutput,"BDToutput/F");
  outTree->Branch("MEoutput",&MEoutput,"MEoutput/F");
  outTree->Branch("ME_p_sig",&p_sig,"ME_p_sig/F");
  outTree->Branch("ME_p_bkg",&p_bkg,"ME_p_bkg/F");
  outTree->Branch("ME_p_err_sig",&p_err_sig,"ME_p_err_sig/F");
  outTree->Branch("ME_p_err_bkg",&p_err_bkg,"ME_p_err_bkg/F");
  outTree->Branch("ME_n_perm_sig",&n_perm_sig,"ME_n_perm_sig/F");
  outTree->Branch("ME_n_perm_bkg",&n_perm_bkg,"ME_n_perm_bkg/F");

 
  outTree->Branch("Weight",&weight,"Weight/F");
  outTree->Branch("Evt_ID",&evtID,"Evt_ID/F");

  BDT = new BDTClassifier("/nfs/dust/cms/user/kelmorab/CMSSW_7_4_15/src/TTH/StandaloneBDT/data/bdtweights_v5");
  MEM = new MEMClassifier();

//   std::cout<<"constructor done"<<std::endl;
  
}

Evaluater::~Evaluater(){
 delete inTree;
 delete outTree;
 delete inTreeFile;
 delete outTreeFile;
//  delete BDTClassifier;
}

void Evaluater::eval(){
  totaltimer->Start();

 long nEvents=inTree->GetEntries();

  char* beginEvtchar = getenv ("BEGINEVENT");
  char* endEvtchar = getenv ("ENDEVENT");

  TString beginEvtstring=TString(beginEvtchar);
  TString endEvtstring=TString(endEvtchar);

  std::cout<<beginEvtstring<<" "<<endEvtstring<<std::endl;

  int beginEvt=0;
  int endEvt=0;

  if(beginEvtstring.IsDigit()){beginEvt=beginEvtstring.Atoi();}
  if(endEvtstring.IsDigit()){endEvt=endEvtstring.Atoi();}
  std::cout<<beginEvt<<" "<<endEvt<<std::endl;

  if(beginEvt==endEvt){std::cout<<"no event number stated, doing all events "<<std::endl; beginEvt=0; endEvt=nEvents-1;}

  for(long ievt=beginEvt;ievt<=endEvt;ievt++){
 //    if(ievt%100==0)std::cout<<"at event "<<ievt<<std::endl;
    std::cout<<"at event "<<ievt<<std::endl;
    inTree->GetEntry(ievt);
    
    std::vector<TLorentzVector> vjetVecs;
    std::vector<double> vjetCSVs;
    TLorentzVector LeptonVec;
    TLorentzVector metVec;
    
    //set up input
//     std::cout<<jetPt[0]<<" "<<jetEta[0]<<" "<<jetPhi[0]<<" "<<jetE[0]<<jetCSV[0]<<std::endl;
    
    vjetVecs=getLVs(nJets, jetPt,jetEta,jetPhi,jetE);
    for(int ijet=0; ijet<nJets;ijet++){
      vjetCSVs.push_back(jetCSV[ijet]);
//       std::cout<<vjetCSVs.back()<<std::endl;
    }
    
    std::vector<TLorentzVector> vLeptonVec;
    std::vector<double> vLeptonCharge;

    LeptonVec=getLV(leptonPt,leptonEta,leptonPhi,leptonE);
    vLeptonVec.push_back(LeptonVec);
    vLeptonCharge.push_back(1.0);
    metVec.SetPtEtaPhiE(metPt,0,metPhi,metPt);
    
//     std::cout<<"eval BDT"<<std::endl;
    Float_t bdtout=BDT->GetBDTOutput(vLeptonVec,vjetVecs,vjetCSVs,vjetVecs,vjetCSVs,metVec);
//     if(bdtout!=oldBDT)std::cout<<bdtout<<" old BDT: "<<oldBDT<<std::endl;
    
    BDToutput=bdtout;
    
//     std::cout<<" BDT Cat "<<BDT->GetCategoryOfLastEvaluation()<<std::endl;
    
//     int ivar=0;
    if(ievt==beginEvt){
      bdtVarMap=BDT->GetVariablesOfLastEvaluation();
      for(std::map<std::string, float >::iterator iter=bdtVarMap.begin(); iter!=bdtVarMap.end(); ++iter) {
//         std::cout << iter->first<<" "<<iter->second << std::endl;
	TString varName=iter->first;
	TString varbrtype=iter->first;
	varbrtype+="/F";
	
	outTree->Branch(varName, &(iter->second),varbrtype);
	}
    }
    else{
      std::map<std::string,float> thisVarMap = BDT->GetVariablesOfLastEvaluation();
      for(std::map<std::string, float >::iterator iter=thisVarMap.begin(); iter!=thisVarMap.end(); ++iter) {
//         std::cout << iter->first<<" "<<iter->second << std::endl;
	bdtVarMap[iter->first]=iter->second;
      }
    }
    
    std::cout<<"BDT done"<<std::endl;
    std::cout<<"BDT result "<<bdtout<<std::endl;
    timer->Start();
    MEMResult thisMEMResult=MEM->GetOutput(vLeptonVec,vLeptonCharge,vjetVecs,vjetCSVs,vjetVecs,vjetCSVs,metVec);
    std::cout<<"mem done"<<std::endl;
    std::cout<<"MEM result "<<thisMEMResult.p<<std::endl;
    std::cout<<"Time for MEM "<<timer->RealTime()<<std::endl;
    
    std::cout<<"BDT vs MEM "<<bdtout<<" "<<thisMEMResult.p<<std::endl;
    
    MEoutput=thisMEMResult.p;

    p_sig=thisMEMResult.p_sig;
    p_bkg=thisMEMResult.p_bkg;
    p_err_sig=thisMEMResult.p_err_sig;
    p_err_bkg=thisMEMResult.p_err_bkg;
    n_perm_sig=thisMEMResult.n_perm_sig;
    n_perm_bkg=thisMEMResult.n_perm_bkg;
    
    outTree->Fill();
//     int pid=getpid();
//     std::cout<<"PID "<<pid<<std::endl;
  }
  	
  outTree->AutoSave();
  outTreeFile->Close();
  
  double totalTime=totaltimer->RealTime();
  std::cout<<"TotalTime "<<totalTime<<std::endl;
} 
  

int main(int argc, char *argv[] ){
  
//   std::cout<<argc<<" "<<argv[1]<<std::endl;
//     std::cout<<"bla"<<std::endl;
// 
//     TString* infileName=new TString(argv[1]);
//     TString* infileName=new TString(argv[1]);
// 
//     std::cout<<"infile name "<<infileName<<std::endl;
  
  char* filename = getenv ("FILENAME");
  char* ofname = getenv ("OUTFILENAME");
  TString infileName=TString(filename);
    TString outfileName=TString(ofname);
    
    Evaluater* myEval=new Evaluater(infileName, outfileName);
    myEval->eval();
    
}
