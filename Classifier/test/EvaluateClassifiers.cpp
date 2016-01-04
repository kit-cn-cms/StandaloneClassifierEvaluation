#include "TTH/CommonClassifier/interface/BDTClassifier.h"
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
#include "TMath.h"


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

  Int_t bdtcategory; 
 
  Float_t jetPt[100]; 
  Float_t jetPhi[100]; 
  Float_t jetEta[100]; 
  Float_t jetE[100]; 
  Float_t jetCSV[100]; 

  Float_t loosejetPt[100]; 
  Float_t loosejetPhi[100]; 
  Float_t loosejetEta[100]; 
  Float_t loosejetE[100]; 
  Float_t loosejetCSV[100]; 
  
  Float_t leptonPt; 
  Float_t leptonPhi; 
  Float_t leptonEta; 
  Float_t leptonE;
  
  Float_t metPt;
  Float_t metPhi;
  
  Int_t nTightMuons;
  Int_t nTightElectrons;
  Float_t muonCharge;
  Float_t electronCharge;
  
  Float_t weight;
  
  Int_t nJets;
  Int_t nlooseJets;

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
  
  inTree->SetBranchAddress("LooseJet_Pt",&loosejetPt);
  inTree->SetBranchAddress("LooseJet_Eta",&loosejetEta);
  inTree->SetBranchAddress("LooseJet_Phi",&loosejetPhi);
  inTree->SetBranchAddress("LooseJet_E",&loosejetE);
  inTree->SetBranchAddress("LooseJet_CSV",&loosejetCSV);

  inTree->SetBranchAddress("N_LooseJets",&nlooseJets);
  inTree->SetBranchAddress("N_Jets",&nJets);
  
  inTree->SetBranchAddress("Evt_Pt_PrimaryLepton",&leptonPt);
  inTree->SetBranchAddress("Evt_Eta_PrimaryLepton",&leptonEta);
  inTree->SetBranchAddress("Evt_Phi_PrimaryLepton",&leptonPhi);
  inTree->SetBranchAddress("Evt_E_PrimaryLepton",&leptonE);
  
  inTree->SetBranchAddress("N_TightElectrons",&nTightElectrons);
  inTree->SetBranchAddress("N_TightMuons",&nTightMuons);
  inTree->SetBranchAddress("Electron_Charge",&electronCharge);
  inTree->SetBranchAddress("Muon_Charge",&muonCharge);
  
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


  outTree->Branch("N_LooseJets",&nlooseJets,"N_LooseJets/I");
  outTree->Branch("N_Jets",&nJets,"N_Jets/I");
  
  outTree->Branch("Jet_Pt",&jetPt,"Jet_Pt[N_Jets]/F");
  outTree->Branch("Jet_Eta",&jetEta,"Jet_Eta[N_Jets]/F");
  outTree->Branch("Jet_Phi",&jetPhi,"Jet_Phi[N_Jets]/F");
  outTree->Branch("Jet_E",&jetE,"Jet_E[N_Jets]/F");
  outTree->Branch("Jet_CSV",&jetCSV,"Jet_CSV[N_Jets]/F");
  
  outTree->Branch("LooseJet_Pt",&loosejetPt,"LooseJet_Pt[N_LooseJets]/F");
  outTree->Branch("LooseJet_Eta",&loosejetEta,"LooseJet_Eta[N_LooseJets]/F");
  outTree->Branch("LooseJet_Phi",&loosejetPhi,"LooseJet_Phi[N_LooseJets]/F");
  outTree->Branch("LooseJet_E",&loosejetE,"LooseJet_E[N_LooseJets]/F");
  outTree->Branch("LooseJet_CSV",&loosejetCSV,"LooseJet_CSV[N_LooseJets]/F");

  
  outTree->Branch("Evt_Pt_PrimaryLepton",&leptonPt,"Evt_Pt_PrimaryLepton/F");
  outTree->Branch("Evt_Eta_PrimaryLepton",&leptonEta,"Evt_Eta_PrimaryLepton/F");
  outTree->Branch("Evt_Phi_PrimaryLepton",&leptonPhi,"Evt_Phi_PrimaryLepton/F");
  outTree->Branch("Evt_E_PrimaryLepton",&leptonE,"Evt_E_PrimaryLepton/F");
  
  outTree->Branch("N_TightElectrons",&nTightElectrons,"N_TightElectrons/I");
  outTree->Branch("N_TightMuons",&nTightMuons,"N_TightMuons/I");
  outTree->Branch("Electron_Charge",&electronCharge,"Electron_Charge/F");
  outTree->Branch("Muon_Charge",&muonCharge,"Muon_Charge/F");
  
  outTree->Branch("Evt_Pt_MET",&metPt,"Evt_Pt_MET/F");
  outTree->Branch("Evt_Phi_MET",&metPhi,"Evt_Phi_MET/F");

  outTree->Branch("BDTcategory",&bdtcategory,"BDTcategory/I");

  outTree->Branch("Weight",&weight,"Weight/F");
  outTree->Branch("Evt_ID",&evtID,"Evt_ID/I");

  BDT = new BDTClassifier("/nfs/dust/cms/user/kelmorab/CMSSW_7_4_15/src/TTH/CommonClassifier/data/bdtweights_v5");
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
  //if(ievt%100==0)std::cout<<"at event "<<ievt<<std::endl;
    std::cout<<"----------------------------------------------------------------------"<<std::endl;
    std::cout<<"at event "<<ievt<<std::endl;
    inTree->GetEntry(ievt);
    
    std::vector<TLorentzVector> vjetVecs;
    std::vector<double> vjetCSVs;
    std::vector<TLorentzVector> vloosejetVecs;
    std::vector<double> vloosejetCSVs;
    TLorentzVector LeptonVec;
    TLorentzVector metVec;
    
    //set up input
//     std::cout<<jetPt[0]<<" "<<jetEta[0]<<" "<<jetPhi[0]<<" "<<jetE[0]<<jetCSV[0]<<std::endl;
    
    vjetVecs=getLVs(nJets, jetPt,jetEta,jetPhi,jetE);
    for(int ijet=0; ijet<nJets;ijet++){
      double thisCSV=TMath::Max( TMath::Min(jetCSV[ijet],1.0), 0.0);
      vjetCSVs.push_back(thisCSV);
//       std::cout<<vjetCSVs.back()<<std::endl;
    }
    
    vloosejetVecs=getLVs(nlooseJets, loosejetPt,loosejetEta,loosejetPhi,loosejetE);
    for(int ijet=0; ijet<nlooseJets;ijet++){
      double thisCSV=TMath::Max( TMath::Min(loosejetCSV[ijet],1.0), 0.0);
      vloosejetCSVs.push_back(thisCSV);
//       std::cout<<vjetCSVs.back()<<std::endl;
    }
    
    std::vector<TLorentzVector> vLeptonVec;
    std::vector<double> vLeptonCharge;

    LeptonVec=getLV(leptonPt,leptonEta,leptonPhi,leptonE);
    vLeptonVec.push_back(LeptonVec);
    
//     std::cout<<nTightElectrons<<" "<<nTightMuons<<" "<<electronCharge<<" "<<muonCharge<<std::endl;
    
    //rewrite this to handel DL case
    if(nTightElectrons==1){
      vLeptonCharge.push_back(electronCharge);
    }
    else if(nTightMuons==1){
      vLeptonCharge.push_back(muonCharge);
    }
    else {
      vLeptonCharge.push_back(1.0);
      std::cout<<"WARNING: nLeptons!=1 -> setting lepton charge to 1.0"<<std::endl;
    }
    metVec.SetPtEtaPhiE(metPt,0,metPhi,metPt);

      
    //print object numbers
    std::cout<<"evtID "<<evtID<<std::endl;
    std::cout<<"Jets pT, eta, phi, m, CSV"<<std::endl;
    for(unsigned int ijet=0; ijet<vjetVecs.size();ijet++){
      std::cout<<vjetVecs.at(ijet).Pt()<<" "<<vjetVecs.at(ijet).Eta()<<" "<<vjetVecs.at(ijet).Phi()<<" "<<vjetVecs.at(ijet).M()<<" "<<vjetCSVs.at(ijet)<<std::endl;
    }
    std::cout<<"Lepton pT, eta, phi, m, charge"<<std::endl;
    for(unsigned int ilep=0; ilep<vLeptonVec.size();ilep++){
      std::cout<<vLeptonVec.at(ilep).Pt()<<" "<<vLeptonVec.at(ilep).Eta()<<" "<<vLeptonVec.at(ilep).Phi()<<" "<<vLeptonVec.at(ilep).M()<<" "<<vLeptonCharge.at(ilep)<<std::endl;
    }
    std::cout<<"MET pT, eta, phi, m"<<std::endl;
    std::cout<<metVec.Pt()<<" "<<metVec.Eta()<<" "<<metVec.Phi()<<" "<<metVec.M()<<std::endl;
    
    
    
//     std::cout<<"eval BDT"<<std::endl;
    Float_t bdtout=BDT->GetBDTOutput(vLeptonVec,vjetVecs,vjetCSVs,vloosejetVecs,vloosejetCSVs,metVec);
//     if(bdtout!=oldBDT)std::cout<<bdtout<<" old BDT: "<<oldBDT<<std::endl;
    
    BDToutput=bdtout;
    TString thisCat=BDT->GetCategoryOfLastEvaluation();
    std::cout<<" BDT Cat "<<BDT->GetCategoryOfLastEvaluation()<<std::endl;

    if(thisCat=="6j2t")bdtcategory=62;
    else if(thisCat=="4j3t")bdtcategory=43;
    else if(thisCat=="5j3t")bdtcategory=53;
    else if(thisCat=="6j3t")bdtcategory=63;
    else if(thisCat=="4j4t")bdtcategory=44;
    else if(thisCat=="5j4t")bdtcategory=54;
    else if(thisCat=="6j4t")bdtcategory=64;
    else bdtcategory=0;
    
    
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
//     timer->Start();
//     MEMResult thisMEMResult=MEM->GetOutput(vLeptonVec,vLeptonCharge,vjetVecs,vjetCSVs,vjetVecs,vjetCSVs,metVec);
//     std::cout<<"mem done"<<std::endl;
//     std::cout<<"MEM result "<<thisMEMResult.p<<std::endl;
//     std::cout<<"Time for MEM "<<timer->RealTime()<<std::endl;
//     
//     std::cout<<"BDT vs MEM "<<bdtout<<" "<<thisMEMResult.p<<std::endl;
//     
//     MEoutput=thisMEMResult.p;
// 
//     p_sig=thisMEMResult.p_sig;
//     p_bkg=thisMEMResult.p_bkg;
//     p_err_sig=thisMEMResult.p_err_sig;
//     p_err_bkg=thisMEMResult.p_err_bkg;
//     n_perm_sig=thisMEMResult.n_perm_sig;
//     n_perm_bkg=thisMEMResult.n_perm_bkg;
    
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
