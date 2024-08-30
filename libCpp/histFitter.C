#ifndef HIST_FITTER
#define HIST_FITTER

#include "TROOT.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"

#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooFFTConvPdf.h"
/// include pdfs
#include "RooCBExGaussShape.h"
#include "RooCMSShape.h"

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <unordered_map>
// #ifdef __CINT__
// #pragma link C++ class std::vector<std::string>+;
// #endif

// Minuit status flags: https://root.cern.ch/doc/master/classROOT_1_1Minuit2_1_1Minuit2Minimizer.html#ab28a14f6c3b1200712e944f986ce63df
// or also https://root-forum.cern.ch/t/meaning-of-values-returned-by-roofitresult-status/16355
// for hesse https://root.cern.ch/doc/v610/classROOT_1_1Minuit2_1_1Minuit2Minimizer.html#afb50781f0567ea8ba364756bdf1d70ad
// Cov matrix quality flags: https://root.cern.ch/doc/master/classROOT_1_1Minuit2_1_1Minuit2Minimizer.html#afc8d76a86981e855014a435a60f3fb84
// for minos https://root.cern.ch/doc/v610/classROOT_1_1Minuit2_1_1Minuit2Minimizer.html#abdff46cdc39c578b941c731ad3b440d1

//using namespace std;
using namespace RooFit;

class tnpFitter {
public:
  tnpFitter( TFile *file, std::string histname, int massbins, float massmin, float massmax  );
  tnpFitter( TH1 *hPass, TH1 *hFail, std::string histname, int massbins, float massmin, float massmax  );
    ~tnpFitter(); //{ if( _work != 0 ) delete _work; }
  void setZLineShapes(TH1 *hZPass, TH1 *hZFail );
  void setTotalBkgShapes(TH1 *hBkgPass, TH1 *hBkgFail );
  void setWorkspace(const std::vector<std::string>&, bool, bool, bool);
  //void setOutputFile(const std::string& fname ) {_fOut = new TFile(fname.c_str(), "recreate"); } 
  void setOutputFile(const std::string& fname );
    //void setOutputFile(const std::string& fname ) {_fname = fname; } 
  void setPlotOutputPath(const std::string& fname) { _outPlotPath = fname;}
  int fits(const std::string& title = "");
  void useMinos(bool minos = true) {_useMinos = minos; }
  void isMC(bool isMC = true) {_isMC = isMC; }
  void textParForCanvas(RooFitResult *resP, RooFitResult *resF, TPad *p, double&, double&);
  void fixSigmaFtoSigmaP(bool fix=true) { _fixSigmaFtoSigmaP= fix; }
  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
  void setZeroBackground(bool zeroBkg = true) {_zeroBackground = zeroBkg; }
  void setPassStrategy(int strategy) {_strategyPassFit = std::clamp(strategy, 0, 2); } 
  void setFailStrategy(int strategy) {_strategyFailFit = std::clamp(strategy, 0, 2); } 
  void setPrintLevel(int level) {_printLevel = std::clamp(level, -1, 9); } 
  void setMaxSignalFractionFail(double max) { _maxSignalFractionFail = max; }
  double getEfficiencyUncertainty(double nP, double nF, double e_nP, double e_nF);
  void updateConstraints(const std::string& key, const std::string& value) { _constraints[key] = value; }
  void setConstantVariable(const std::string& name, const double& val, const bool& removeRange);
  RooFitResult* manageFit(bool, int, std::string*, double*);
    
private:
  RooWorkspace *_work;
  std::string _histname_base = "";
  TFile *_fOut;
  std::string _fname = "";
  double _nTotP, _nTotF;
  bool _useMinos = false;
  bool _isMC = false;
  bool _zeroBackground = false;  
  bool _fixSigmaFtoSigmaP = false;
  double _xFitMin,_xFitMax;
  int _strategyPassFit = 1;
  int _strategyFailFit = 1;
  int _printLevel = 3;
  double _maxSignalFractionFail = -1; // not used by default if negative
  std::unordered_map<std::string, std::string> _constraints = {};
  int _nFitBins = -1;
  bool _hasShape_bkgFailMC = false;
  std::string _outPlotPath = "";
    
};

tnpFitter::tnpFitter(TFile *filein, std::string histname, int massbins, float massmin, float massmax ) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;  

  std::string namePass = TString::Format("%s_Pass",histname.c_str()).Data();
  std::string nameFail = TString::Format("%s_Fail",histname.c_str()).Data();
  TH1 *hPass = (TH1*) filein->Get(namePass.c_str());
  TH1 *hFail = (TH1*) filein->Get(nameFail.c_str());
  if (hPass == nullptr) {
      std::cout << "Error reading " << namePass << " from " << filein->GetName() << std::endl;
      exit(EXIT_FAILURE);
  }
  if (hFail == nullptr) {
      std::cout << "Error reading " << nameFail << " from " << filein->GetName() << std::endl;
      exit(EXIT_FAILURE);
  }
  
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ ) {
      if(  hPass->GetXaxis()->GetBinCenter(ib) <= massmin || hPass->GetXaxis()->GetBinCenter(ib) >= massmax ) {
          hPass->SetBinContent(ib,0);
          hFail->SetBinContent(ib,0);
      }
      // protection for Chi2
      if (hPass->GetBinError(ib) <= 0.0) hPass->SetBinError(ib, 1.0);
      if (hFail->GetBinError(ib) <= 0.0) hFail->SetBinError(ib, 1.0);
  }
      
  _work = new RooWorkspace("w") ;
  //_work->factory("x[50,130]");
  _work->factory(TString::Format("x[%f,%f]",massmin, massmax));

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;
  _xFitMin = massmin;
  _xFitMax = massmax;
  _nFitBins = massbins;
}

tnpFitter::tnpFitter(TH1 *hPass, TH1 *hFail, std::string histname, int massbins, float massmin, float massmax ) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  _histname_base = histname;
  
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ ) {
      if (hPass->GetXaxis()->GetBinCenter(ib) < massmin or hPass->GetXaxis()->GetBinCenter(ib) > massmax) {
          hPass->SetBinContent(ib,0);
          hFail->SetBinContent(ib,0);
      }
      // protection for Chi2
      if (hPass->GetBinError(ib) <= 0.0) hPass->SetBinError(ib, 1.0);
      if (hFail->GetBinError(ib) <= 0.0) hFail->SetBinError(ib, 1.0);
  }      
 
  _work = new RooWorkspace("w") ;
  //_work->factory("x[50,130]");
  _work->factory(TString::Format("x[%f,%f]",massmin, massmax));

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;
  _xFitMin = massmin;
  _xFitMax = massmax;
  _nFitBins = massbins;
  
}

tnpFitter::~tnpFitter() {
    const char* fname = _fOut->GetName();
    if( _work != 0 )
        delete _work;
    if (_fOut and _fOut->IsOpen()) {
        // std::cout << ">>> Closing file " << fname << std::endl;
        _fOut->Close();
    }

    // check goodness of file inside this job     
    //std::cout << "Inside destructor: check goodness of file " << fname << std::endl;
    bool isGood = true;
    TFile* fcheck = TFile::Open(fname, "READ");
    if (not fcheck or fcheck->IsZombie()) {
        isGood = false;
    } else {
        if (fcheck->GetSize() < 512) isGood = false; // set limit at 0.5 kB, file is actually larger
        //if (fcheck->GetSize() < 500000) isGood = false; // set limit at 500 kB, file is actually larger
        else if (fcheck->TestBit(TFile::kRecovered)) isGood = false;
        fcheck->Close();
    }
  
    if (isGood) {
        //std::cout << ">>>> File looks good" << std::endl;
    } else {
        std::cout << "#### File is bad or non existing. Will be deleted if existing" << std::endl;
        if (not gSystem->AccessPathName(fname)) {
            // file exists, let's delete it   
            gSystem->Unlink(fname); // this works also for non-Unix systems, just in case
        }
    }
    
}

//// Let's open this file directly before using it, to avoid root sorceries with memory and object ownership
void tnpFitter::setOutputFile(const std::string& fname) {
    _fOut = new TFile(fname.c_str(), "recreate");
    if (!_fOut || _fOut->IsZombie()) {
        std::cout << "Error opening file " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
}

void tnpFitter::setConstantVariable(const std::string& name, const double& val = 0.0, const bool& removeRange = false) {
    RooRealVar* tmp = _work->var(name.c_str());
    if (tmp != nullptr) {
        if (removeRange) tmp->removeRange();
        tmp->setVal(val);
        tmp->setConstant();
    }
}


void tnpFitter::setZLineShapes(TH1 *hZPass, TH1 *hZFail ) {
  RooDataHist rooPass("hGenZPass","hGenZPass",*_work->var("x"),hZPass);
  RooDataHist rooFail("hGenZFail","hGenZFail",*_work->var("x"),hZFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;  
}

void tnpFitter::setTotalBkgShapes(TH1 *hBkgPass, TH1 *hBkgFail ) {
  RooDataHist rooPass("hTotBkgPass","hTotBkgPass",*_work->var("x"),hBkgPass);
  RooDataHist rooFail("hTotBkgFail","hTotBkgFail",*_work->var("x"),hBkgFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;  
}

void tnpFitter::setWorkspace(const std::vector<std::string>& workspace, bool isMCfit = false, bool analyticPhysicsShape = false, bool modelFSR = false) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }

  if (not analyticPhysicsShape) {
      _work->factory("HistPdf::sigPhysPass(x,hGenZPass,3)");
      _work->factory("HistPdf::sigPhysFail(x,hGenZFail,3)");
  }
  // this x variable should only be needed for the convolution, since the actual binning comes from the histograms
  // increase number of bins and also the range so to span the whole range where the pdfs is larger than 0 (maybe the range is less important here)
  // see also https://root-forum.cern.ch/t/bad-fit-at-boundaries-for-convoluted-roohistpdf/21980/9
  _work->var("x")->setBins(2000, "cache"); // sometimes 10k works, but in newer root version 10k is the maximum including the buffer apparently
  // _work->var("x")->setMin("cache", 50.0); 
  // _work->var("x")->setMax("cache", 130.0); 
  _work->factory(TString::Format("nSigP[%f,0.5,%f]",_nTotP*0.9,_nTotP*1.5));
  RooFFTConvPdf* convPass = (RooFFTConvPdf*) _work->factory("FCONV::sigPass(x, sigPhysPass , sigResPass)");
  convPass->setBufferFraction(0.5);
  
  if (_zeroBackground) {
      // to implement properly
      _work->factory("nBkgP[0]");
      std::cout << "Setting background to zero for pass pdf" << std::endl;
  } else {
      _work->factory(TString::Format("nBkgP[%f,0.5,%f]",_nTotP*0.1,_nTotP*1.5));
  }
  _work->factory("SUM::pdfPass(nSigP*sigPass,nBkgP*bkgPass)");
  
  if (modelFSR) {

      if (isMCfit) {
          _work->factory(TString::Format("nSigF[%f,%f,%f]",_nTotF*0.9,_nTotF*0.85,_nTotF*1.5));
          if (_zeroBackground) {
              // to implement properly
              _work->factory("nBkgF[0]");
              std::cout << "Setting background to zero for fail pdf" << std::endl;
          } else {
              _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*0.15));
          }
      } else{ 
          _work->factory(TString::Format("nSigF[%f,0.5,%f]",_nTotF*0.9,_nTotF*1.5));
          _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*1.5));
      }

      RooFFTConvPdf* convFail = (RooFFTConvPdf*) _work->factory("FCONV::sigMainFail(x, sigPhysFail , sigResFail)");
      convFail->setBufferFraction(0.5);
      _work->factory("SUM::sigFail(fracMainF[0.95,0.8,1.0]*sigMainFail, sigFsrFail)");
      _work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
          
  } else {

      if (isMCfit) {
          _work->factory(TString::Format("nSigF[%f,%f,%f]",_nTotF*0.9,_nTotF*0.85,_nTotF*1.5));
          if (_zeroBackground) {
              // to implement properly
              _work->factory("nBkgF[0]");
              std::cout << "Setting background to zero for fail pdf" << std::endl;
          } else {
              _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*0.15));
          }
      } else{ 
          if (_work->var("maxFracSigF") != nullptr) {
              // std::cout << "Test setting signal fraction" << std::endl;
              double maxFracSigF = _work->var("maxFracSigF")->getVal(); 
              double minFracBkgF = 1.0 - maxFracSigF;
              double halfFracSigF = maxFracSigF/2.0;
              _work->factory(TString::Format("nSigF[%f,0.5,%f]", _nTotF*halfFracSigF, _nTotF*maxFracSigF));
              _work->factory(TString::Format("nBkgF[%f,%f,%f]", _nTotF*(minFracBkgF+halfFracSigF), _nTotF*minFracBkgF,_nTotF*1.5));          
          } else {
              _work->factory(TString::Format("nSigF[%f,0.5,%f]",_nTotF*0.9,_nTotF*1.5));
              _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*1.5));          
          }
      }
      RooFFTConvPdf* convFail = (RooFFTConvPdf*) _work->factory("FCONV::sigFail(x, sigPhysFail , sigResFail)");
      convFail->setBufferFraction(0.5);
      _work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
      
  }

  if (_work->pdf("bkgFailBackup") != nullptr)
      _work->factory("SUM::pdfFailBackup(nSigF*sigFail,nBkgF*bkgFailBackup)");
  if (_isMC and _work->pdf("bkgFailMC") != nullptr) {
      _work->factory("SUM::pdfFailMC(nSigF*sigFail,nBkgF*bkgFailMC)");
      _hasShape_bkgFailMC = true;
  }
  //_work->Print(); // FIXME: might want to comment this one to avoid unnecessary output text

}

RooFitResult* tnpFitter::manageFit(bool isPass, int attempt = 0, std::string* lastNamePDF = nullptr, double* chi2value = nullptr) {
    
    std::string pdfName = "pdfPass";
    std::string hName = "hPass";
    std::string constrainName = "constrainP";
    std::string sigPar = "nSigP";
    std::string bkgPar = "nBkgP";

    if (not isPass) {
        if (_isMC and _hasShape_bkgFailMC) {
            pdfName = (attempt == 2) ? "pdfFailBackup" : "pdfFailMC";
        } else {
            pdfName = (attempt == 2) ? "pdfFailBackup" : "pdfFail";
        }
        hName = "hFail";
        constrainName = "constrainF";
        sigPar = "nSigF";
        bkgPar = "nBkgF";
    }

    if (lastNamePDF != nullptr) *lastNamePDF = pdfName;

    // should check the parameters of the actual pdf being used rather than assuming that there are no constraints when attempt == 2
    const RooArgSet* constraint = (attempt == 2) ? nullptr : _work->set(constrainName.c_str());
    if (_isMC and _hasShape_bkgFailMC) constraint = nullptr;
    // TODO: in case minos is used the uncertainties are asymmetric and one should get them accordingly
    
    RooAbsPdf *pdf = _work->pdf(pdfName.c_str());
    //RooAbsData* dh = _work->data(hName.c_str());
    RooAbsData* dh =  _work->data(hName.c_str());
    RooFitResult* res = pdf->fitTo(*dh,
                                   // Minos(_work->argSet(sigPar.c_str())),
                                   _useMinos ? Minos(_work->argSet(sigPar.c_str())) : Minos(kFALSE),
                                   (_useMinos or not _isMC) ? SumW2Error(kFALSE) : SumW2Error(kTRUE), // try always false
                                   //SumW2Error(kFALSE), // default is false, but needs it explicitly for MC otherwise roofit complains
                                   Save(),
                                   Range("fitMassRange"),
                                   Minimizer("Minuit2"),
                                   Strategy(isPass ? _strategyPassFit : _strategyFailFit),
                                   PrintLevel(_printLevel),
                                   (constraint != nullptr) ? ExternalConstraints(*constraint) : RooCmdArg::none() );

    RooChi2Var chi2("chi2", "chi2 var", *pdf, *((RooDataHist*) dh), Range(_xFitMin,_xFitMax));
    *chi2value = chi2.getVal();
    int ndof = _nFitBins - res->floatParsFinal().getSize();
    double chi2sigma = std::sqrt(2. * ndof);
    bool goodChi2 = std::fabs(*chi2value - (double) ndof) < (10.0 * chi2sigma); 
    // std::cout << pdfName << " --> Chi2 / ndof = " << chi2.getVal() << " / " << ndof << "   good chi2 = " << goodChi2 << std::endl;
    // std::cout << pdfName << " --> status = " << res->status() << std::endl;
    // std::cout << pdfName << " --> cov. quality = " << res->covQual() << std::endl;

    if (attempt > 0) return res;

    if (_isMC) {
        if (goodChi2 and (res->status() == 0 or res->status() == 1)) return res;
    } else {
        if (goodChi2 and res->covQual() == 3 and (res->status() == 0 or res->status() == 1)) return res;
    }
    //std::cout << "Failed fit for " << pdfName << ": trying again ..." << std::endl;
    // if status != 0 try something, like checking background and if it is too small fit only signal
    // or just refit with ranges of parameters frm previous fit, but this might not work
    double nBkg = _work->var(bkgPar.c_str())->getVal();
    double nSig = _work->var(sigPar.c_str())->getVal();
    RooAbsPdf *bkgpdf = _work->pdf(isPass ? "bkgPass" : "bkgFail");
    
    if (nBkg < 0.005 * nSig) {
        // a bit hardcoded to get background parameters, should probably make sure to have these stored somewhere in the class
        //std::cout << "Refitting with no background: attempt " << attempt+1 << std::endl;
        setConstantVariable(bkgPar.c_str(), 0.0, true);
        std::vector<std::string> bkgParNamesP = {"acmsP", "betaP", "gammaP", "expalphaP"};
        std::vector<std::string> bkgParNamesF = {"acmsF", "betaF", "gammaF", "expalphaF"};
        std::vector<std::string>& bkgParNames = isPass ? bkgParNamesP : bkgParNamesF; 
        for (UInt_t i = 0; i < bkgParNames.size(); i++) {
            if (_work->var(bkgParNames[i].c_str()) == nullptr) continue;
            setConstantVariable(bkgParNames[i], _work->var(bkgParNames[i].c_str())->getVal());
        }
        return manageFit(isPass, 1, lastNamePDF, chi2value);
    } else {
        if (isPass) return res;
        else        return manageFit(isPass, 2, lastNamePDF, chi2value);
    }
    
}


int tnpFitter::fits(const std::string& title) {

  // std::cout << " this is the title : " << title << std::endl;

  
  // RooAbsPdf *pdfPass = _work->pdf("pdfPass");
  // RooAbsPdf *pdfFail = _work->pdf("pdfFail");

  /// FC: seems to be better to change the actual range than using a fitRange in the fit itself (???)
  /// FC: I don't know why but the integral is done over the full range in the fit not on the reduced range
  _work->var("x")->setRange(_xFitMin,_xFitMax);
  _work->var("x")->setRange("fitMassRange",_xFitMin,_xFitMax);

  // TODO: check that all parameters exists
  for (auto i = _constraints.begin(); i != _constraints.end(); i++) {
      // std::cout << i->first << " -> " << i->second << std::endl;
      _work->defineSet((i->first).c_str(), (i->second).c_str());
  }

  std::string lastNamePassPDF = "";
  double chi2valuePass = 0.0;
  //std::cout << "Fit for passing probes" << std::endl;
  RooFitResult* resPass = manageFit(true, 0, &lastNamePassPDF, &chi2valuePass);
  // std::cout << "Last used pass pdf  -> " << lastNameFailPDF << std::endl;
  
  //std::cout << "Fit for failing probes" << std::endl;
  std::string lastNameFailPDF = "";
  double chi2valueFail = 0.0;
  RooFitResult* resFail = manageFit(false, 0, &lastNameFailPDF, &chi2valueFail);
  // std::cout << "Last used fail pdf  -> " << lastNameFailPDF << std::endl;

  std::string bkgNamePass = "bkgPass";
  std::string bkgNameFail = "bkgFail";
  if (lastNamePassPDF.find("Backup") != std::string::npos) bkgNamePass += "Backup";
  if (lastNameFailPDF.find("Backup") != std::string::npos) {
      bkgNameFail += "Backup";
  } else {
      if (_isMC and _hasShape_bkgFailMC) bkgNameFail = "bkgFailMC"; 
  }
      
  RooPlot *pPass = _work->var("x")->frame(_xFitMin,_xFitMax); // always plot 50 - 130
  RooPlot *pFail = _work->var("x")->frame(_xFitMin,_xFitMax);
  pPass->SetTitle("passing probe");
  pFail->SetTitle("failing probe");

  _work->data("hPass") ->plotOn( pPass, Name("data_pass") );
  _work->pdf(lastNamePassPDF.c_str())->plotOn( pPass, LineColor(kRed), Name("model_pass") );
  _work->pdf(lastNamePassPDF.c_str())->plotOn( pPass, Components(bkgNamePass.c_str()),LineColor(kBlue),LineStyle(kDashed), Name("bkg_pass"));
  _work->data("hPass") ->plotOn( pPass );
  
  _work->data("hFail") ->plotOn( pFail, Name("data_fail") );
  _work->pdf(lastNameFailPDF.c_str())->plotOn( pFail, LineColor(kRed), Name("model_fail") );
  _work->pdf(lastNameFailPDF.c_str())->plotOn( pFail, Components(bkgNameFail.c_str()),LineColor(kBlue),LineStyle(kDashed),Name("bkg_fail"));
  _work->data("hFail") ->plotOn( pFail );

  std::string canvasName = _histname_base + "_Canv"; // TString::Format("%s_Canv",_histname_base.c_str()); 
  TCanvas * c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1150, 500);
  c->Divide(3,1);
  TPad *padText = (TPad*)c->GetPad(1);
  textParForCanvas( resPass,resFail, padText, chi2valuePass, chi2valueFail);
  c->cd(2);
  pPass->Draw();
  c->cd(3);
  pFail->Draw();

  c->cd(0);
  c->SaveAs(TString::Format("%s%s.png", _outPlotPath.c_str(), canvasName.c_str())); // not sure .pdf would save all 3 pads, maybe not that simply
  
  // _fOut = new TFile(_fname.c_str(), "recreate");
  // if (!_fOut || _fOut->IsZombie()) {
  //     std::cout << "Error opening file " << _fname << std::endl;
  //     exit(EXIT_FAILURE);
  // }
  _fOut->cd();
  c->Write(canvasName.c_str(), TObject::kOverwrite);
  // pPass->Write(TString::Format("%s_rooplotP", _histname_base.c_str()), TObject::kOverwrite);
  // pFail->Write(TString::Format("%s_rooplotF", _histname_base.c_str()), TObject::kOverwrite);
  resPass->Write(TString::Format("%s_resP",_histname_base.c_str()), TObject::kOverwrite);
  resFail->Write(TString::Format("%s_resF",_histname_base.c_str()), TObject::kOverwrite);
  //_fOut->Close(); // closed in the destructor
  
  return 1;

}



double tnpFitter::getEfficiencyUncertainty(double nP, double nF, double e_nP, double e_nF) {

    double nTot = nP + nF; 
    return 1./(nTot*nTot) * std::sqrt( nP*nP* e_nF*e_nF + nF*nF * e_nP*e_nP );
    
}

/////// Stupid parameter dumper /////////
void tnpFitter::textParForCanvas(RooFitResult *resP, RooFitResult *resF,TPad *p, double& chi2valuePass, double& chi2valueFail) {

  double eff = -1;
  double e_eff = 0;

  RooRealVar *nSigP = _work->var("nSigP");
  RooRealVar *nSigF = _work->var("nSigF");
  
  double nP   = nSigP->getVal();
  double e_nP = nSigP->getError();
  double nF   = nSigF->getVal();
  double e_nF = nSigF->getError();
  double nTot = nP+nF;
  eff = nP / (nP + nF);
  e_eff = getEfficiencyUncertainty(nP, nF, e_nP, e_nF);  // this is linear error propagation assuming uncorrelated nP and nF, but might not be correct when efficiency is close to 1

  // nP and nF should always have uncertainties equal to at least sqrt(n) in real data
  double e_eff_corr = e_eff;
  if (not _isMC) {
      double e_nP_corr = std::max(e_nP, std::sqrt(nP));
      double e_nF_corr = std::max(e_nF, std::sqrt(nF));
      // std::cout << "Corrected stat uncertainties on nP and nF --> " << e_nP_corr << ", " << e_nF_corr << std::endl;
      e_eff_corr = getEfficiencyUncertainty(nP, nF, e_nP_corr, e_nF_corr); 
  }

  
  TPaveText *text1 = new TPaveText(0,0.76,1,1);
  text1->SetFillColor(0);
  text1->SetBorderSize(0);
  text1->SetTextAlign(12);

  // better to just print the status at the end, it is what really matters and it is less lines to print
  // for (UInt_t i = 0 ; i < resP->numStatusHistory(); i++) {
  //     text1->AddText(TString::Format("%s status: pass %d, fail %d", resP->statusLabelHistory(i), resP->statusCodeHistory(i), resF->statusCodeHistory(i)));
  // }
  text1->AddText(TString::Format("fit status:  pass %d, fail %d",resP->status(),resF->status()));
  text1->AddText(TString::Format("cov quality: pass %d, fail %d",resP->covQual(),resF->covQual()));
  int ndofP = _nFitBins - resP->floatParsFinal().getSize();
  int ndofF = _nFitBins - resF->floatParsFinal().getSize();
  double chi2probPass = 100.0 * TMath::Prob(chi2valuePass, ndofP);
  double chi2probFail = 100.0 * TMath::Prob(chi2valueFail, ndofF);
  text1->AddText(TString::Format("#Chi^{2} (prob): P %.1f/%d (%.1f%%), F %.1f/%d (%.1f%%)", chi2valuePass, ndofP, chi2probPass, chi2valueFail, ndofF, chi2probFail));
  //text1->SetTextFont(62);
  if (not _isMC and (e_eff_corr > e_eff) ) {
      text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f (%1.4f)",eff,e_eff, e_eff_corr));
  } else {
      text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f",eff,e_eff));
  }
  
  //  text->SetTextSize(0.06);

//  text->AddText("* Passing parameters");
  TPaveText *text = new TPaveText(0,0,1,0.76);
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  std::string bkgWarning = "";
  // std::cout << "nBkgP " << _work->var("nBkgP")->getVal() << " +/- " << _work->var("nBkgP")->getError() << std::endl;
  // std::cout << "nBkgF " << _work->var("nBkgF")->getVal() << " +/- " << _work->var("nBkgF")->getError() << std::endl;
  if (_work->var("nBkgP")->getVal() <= 0.0) bkgWarning += "  nBkgP=0";
  if (_work->var("nBkgF")->getVal() <= 0.0) bkgWarning += "  nBkgF=0";
 
  text->AddText(TString::Format("    --- parameters %s", bkgWarning.c_str()) );

  RooArgList listParFinalP = resP->floatParsFinal();
  for( int ip = 0; ip < listParFinalP.getSize(); ip++ ) {
    TString vName = listParFinalP[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }

//  text->AddText("* Failing parameters");
  RooArgList listParFinalF = resF->floatParsFinal();
  for( int ip = 0; ip < listParFinalF.getSize(); ip++ ) {
    TString vName = listParFinalF[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }

  p->cd();
  text1->Draw();
  text->Draw();
}


#endif
