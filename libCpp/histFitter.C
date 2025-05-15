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
#include "RooFFTConvPdf.h"
/// include pdfs
#include "RooGaussian.h"
#include "RooCBExGaussShape.h"
#include "RooCMSShape.h"
#include "RooParamHistFunc.h"

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
    void setBarlowBeestonBkgPdf(bool isPass);
    void evalChi2(const std::string& pdfName, const std::string& hName, const int nFloatPars, bool isPass);
    void setWorkspace(const std::vector<std::string>&, bool, bool, bool);
    //void setOutputFile(const std::string& fname ) {_fOut = new TFile(fname.c_str(), "recreate"); } 
    void setOutputFile(const std::string& fname );
    //void setOutputFile(const std::string& fname ) {_fname = fname; } 
    void setPlotOutputPath(const std::string& fname) { _outPlotPath = fname;}
    int fits(const std::string& title = "");
    void useMinos(bool minos = true) {_useMinos = minos; }
    void isMC(bool isMC = true) {_isMC = isMC; }
    void textParForCanvas(TPad *p, RooFitResult *resP, RooFitResult *resF, double&, double&);
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
    int _nFitBins = -1;
    double _chi2P, _chi2F;
    int _ndofP, _ndofF;
    double _maxSignalFractionFail = -1; // not used by default if negative
    std::unordered_map<std::string, std::string> _constraints = {};
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
    for (int ib=0; ib<=hPass->GetXaxis()->GetNbins()+1; ib++) {
        if (hPass->GetXaxis()->GetBinCenter(ib) <= massmin || hPass->GetXaxis()->GetBinCenter(ib) >= massmax) {
            hPass->SetBinContent(ib,0);
            hFail->SetBinContent(ib,0);
        }
        // protection for Chi2
        if (hPass->GetBinError(ib) <= 0.0) hPass->SetBinError(ib,1.0);
        if (hFail->GetBinError(ib) <= 0.0) hFail->SetBinError(ib,1.0);
    }
      
    _work = new RooWorkspace("w");
    _work->factory(TString::Format("x[%f,%f]", massmin, massmax));

    RooDataHist rooPass("hPass", "hPass", *_work->var("x"), hPass);
    RooDataHist rooFail("hFail", "hFail", *_work->var("x"), hFail);
    _work->import(rooPass);
    _work->import(rooFail);
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
    for (int ib=0; ib<=hPass->GetXaxis()->GetNbins()+1; ib++) {
        if (hPass->GetXaxis()->GetBinCenter(ib) <= massmin || hPass->GetXaxis()->GetBinCenter(ib) >= massmax) {
            hPass->SetBinContent(ib,0);
            hFail->SetBinContent(ib,0);
        }
        // protection for Chi2
        if (hPass->GetBinError(ib) <= 0.0) hPass->SetBinError(ib,1.0);
        if (hFail->GetBinError(ib) <= 0.0) hFail->SetBinError(ib,1.0);
    }      

    _work = new RooWorkspace("w");
    _work->factory(TString::Format("x[%f,%f]", massmin, massmax));

    RooDataHist rooPass("hPass", "hPass", *_work->var("x"), hPass);
    RooDataHist rooFail("hFail", "hFail", *_work->var("x"), hFail);
    _work->import(rooPass);
    _work->import(rooFail);
    _xFitMin = massmin;
    _xFitMax = massmax;
    _nFitBins = massbins;
  
}

tnpFitter::~tnpFitter() {
    const char* fname = _fOut->GetName();
    if (_work != 0)
        delete _work;
    if (_fOut && _fOut->IsOpen()) {
        _fOut->Close();
    }

    // check goodness of file inside this job     
    bool isGood = true;
    TFile* fcheck = TFile::Open(fname, "READ");
    if (!fcheck || fcheck->IsZombie()) {
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
        if (!gSystem->AccessPathName(fname)) {
            gSystem->Unlink(fname); // this works also for non-Unix systems, just in case
        }
    }
    
}

void tnpFitter::setOutputFile(const std::string& fname) {
    _fOut = new TFile(fname.c_str(), "RECREATE");
    if (!_fOut || _fOut->IsZombie()) {
        std::cout << "Error opening file " << fname << std::endl;
        exit(EXIT_FAILURE);
    }
}

void tnpFitter::setConstantVariable(const std::string& name, const double& val=0.0, const bool& removeRange=false) {
    RooRealVar* tmp = _work->var(name.c_str());
    if (tmp != nullptr) {
        if (removeRange) tmp->removeRange();
        tmp->setVal(val);
        tmp->setConstant();
    }
}

void tnpFitter::setZLineShapes(TH1 *hZPass, TH1 *hZFail) {
    RooDataHist rooPass("hZPass", "hZPass", *_work->var("x"), hZPass);
    RooDataHist rooFail("hZFail", "hZFail", *_work->var("x"), hZFail);
    _work->import(rooPass);
    _work->import(rooFail);  
}

void tnpFitter::setTotalBkgShapes(TH1 *hBkgPass, TH1 *hBkgFail) {
    RooDataHist rooPass("hBkgPass", "hBkgPass", *_work->var("x"), hBkgPass);
    RooDataHist rooFail("hBkgFail", "hBkgFail", *_work->var("x"), hBkgFail);
    _work->import(rooPass);
    _work->import(rooFail);  
}

void tnpFitter::setBarlowBeestonBkgPdf(bool isPass=false) {
    std::string hName = isPass ? "hBkgPass" : "hBkgFail";
    std::string pdfName = isPass ? "bkgPass" : "bkgFail";  
    std::string paramHistName = isPass ? "paramHistP" : "paramHistF";

    _work->factory("one[1]");
    _work->factory(TString::Format("RooParamHistFunc::%s(%s,x)", paramHistName.c_str(), hName.c_str()));
    _work->factory(TString::Format("RooRealSumPdf::%s(%s,one)", pdfName.c_str(), paramHistName.c_str()));
}

void tnpFitter::evalChi2(const std::string& pdfName, const std::string& hName, const int nFloatPars, bool isPass=true) {
    RooAbsPdf* pdf = _work->pdf(pdfName.c_str());
    RooAbsData* dh =  _work->data(hName.c_str());
    RooRealVar* x = _work->var("x");
    double sumLL = 0.0, maxLL = 0.0;
    int usedBins = 0;

    std::string parHist = isPass ? "paramHistP" : "paramHistF";
    std::string constrainCategory = isPass ? "constrainP" : "constrainF";
    const RooArgSet* constraint = _work->set(constrainCategory.c_str());

    double binVolume = (_xFitMax-_xFitMin)/_nFitBins;

    for (int ib=0; ib<dh->numEntries(); ib++) {
        x->setVal(_work->var("x")->getBinning().binCenter(ib));
        const RooArgSet* xSet = dh->get(ib);
        double weight = dh->weight();
        double pdfval = pdf->getVal(RooArgSet(*_work->var("x")));
        double mu = std::max(dh->sumEntries()*pdfval*binVolume, 0.1);

        if (weight > 0) {
            sumLL += 2*(weight*TMath::Log(mu) - mu);
            maxLL += 2*(weight*TMath::Log(weight) - weight);
            usedBins++;

            if (_work->function(parHist.c_str()) != nullptr) {

                RooRealVar* nFittedBkg = _work->var((isPass) ? "nBkgP" : "nBkgF");
                RooAbsData* dhBkg = _work->data((isPass) ? "hBkgPass" : "hBkgFail");
                RooParamHistFunc* paramHist = (RooParamHistFunc*)_work->function(parHist.c_str());

                const RooArgSet* xSet = dhBkg->get(ib);
                double weight_bkg = dhBkg->weight() * (nFittedBkg->getVal()/dhBkg->sumEntries());
                
                std::string gammaName = parHist + "_gamma_bin_" + std::to_string(ib);
                RooRealVar* gammaVar = dynamic_cast<RooRealVar*>(*_work->function(parHist.c_str())->servers().findByName(gammaName.c_str()));
                double gamma = gammaVar->getVal();

                sumLL += 2*(weight_bkg*TMath::Log(weight_bkg*gamma) - weight_bkg*gamma);
                maxLL += 2*(weight_bkg*TMath::Log(weight_bkg) - weight_bkg);
            }
        }
    }

    if (constraint != nullptr) {
        for (auto constrPdf=constraint->begin(); constrPdf != constraint->end(); ++constrPdf) {
            RooGaussian* constraintPdf = dynamic_cast<RooGaussian*>(*constrPdf);
            if (constraintPdf && constraintPdf->InheritsFrom("RooGaussian")) {
                double mean = constraintPdf->getMean().getVal();
                double sigma = constraintPdf->getSigma().getVal();
                std::string fitParName = std::string(constraintPdf->GetName()).erase(0, 11); // remove "constraintP_" or "constraintF_"
                double fitVal = _work->var(fitParName)->getVal();

                sumLL -= (fitVal - mean)*(fitVal - mean)/(sigma*sigma);
            }
        }
    }

    double& chi2 = (isPass) ? _chi2P : _chi2F;
    int& ndof    = (isPass) ? _ndofP : _ndofF;

    chi2 = maxLL - sumLL;
    ndof = usedBins - nFloatPars;
           
}

void tnpFitter::setWorkspace(const std::vector<std::string>& workspace, bool isMCfit=false, bool analyticPhysicsShape=false, bool modelFSR=false) {

    for (unsigned icom=0; icom<workspace.size(); ++icom) {
        _work->factory(workspace[icom].c_str());
    }

    if (!analyticPhysicsShape) {
        _work->factory("HistPdf::sigPhysPass(x,hZPass,3)");
        _work->factory("HistPdf::sigPhysFail(x,hZFail,3)");
    }
    // this x variable should only be needed for the convolution, since the actual binning comes from the histograms
    // increase number of bins and also the range so to span the whole range where the pdfs is larger than 0 (maybe the range is less important here)
    // see also https://root-forum.cern.ch/t/bad-fit-at-boundaries-for-convoluted-roohistpdf/21980/9
    _work->var("x")->setBins(2000, "cache"); // sometimes 10k works, but in newer root version 10k is the maximum including the buffer apparently
    // _work->var("x")->setMin("cache", 50.0); 
    // _work->var("x")->setMax("cache", 130.0);

    _work->factory(TString::Format("nSigP[%f,0.5,%f]", _nTotP*0.9, _nTotP*1.5));
    RooFFTConvPdf* convPass = (RooFFTConvPdf*)_work->factory("FCONV::sigPass(x,sigPhysPass,sigResPass)");
    convPass->setBufferFraction(0.5);

    if (_zeroBackground) {
        _work->factory("nBkgP[0]");
        std::cout << "Setting background to zero for pass pdf" << std::endl;
    } else {
        _work->factory(TString::Format("nBkgP[%f,0.5,%f]", _nTotP*0.1, _nTotP*1.5));
    }

    _work->factory("SUM::pdfPass(nSigP*sigPass, nBkgP*bkgPass)");

    if (modelFSR) {

        if (isMCfit) {
            _work->factory(TString::Format("nSigF[%f,%f,%f]",_nTotF*0.9,_nTotF*0.85,_nTotF*1.5));
            if (_zeroBackground) {
                // to implement properly
                _work->factory("nBkgF[0]");
                std::cout << "Setting background to zero for fail pdf" << std::endl;
            } else {
                _work->factory(TString::Format("nBkgF[%f,0.5,%f]", _nTotF*0.1, _nTotF*0.15));
            }
        } else { 
            _work->factory(TString::Format("nSigF[%f,0.5,%f]", _nTotF*0.9, _nTotF*1.5));
            _work->factory(TString::Format("nBkgF[%f,0.5,%f]", _nTotF*0.1, _nTotF*1.5));
        }

        RooFFTConvPdf* convFail = (RooFFTConvPdf*) _work->factory("FCONV::sigMainFail(x, sigPhysFail, sigResFail)");
        convFail->setBufferFraction(0.5);
        _work->factory("SUM::sigFail(fracMainF[0.95,0.8,1.0]*sigMainFail, sigFsrFail)");
        _work->factory("SUM::pdfFail(nSigF*sigFail, nBkgF*bkgFail)");
            
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
        } else { 
            if (_work->var("maxFracSigF") != nullptr) {
                // std::cout << "Test setting signal fraction" << std::endl;
                double maxFracSigF = _work->var("maxFracSigF")->getVal(); 
                double minFracBkgF = 1.0 - maxFracSigF;
                double halfFracSigF = maxFracSigF/2.0;
                _work->factory(TString::Format("nSigF[%f,0.5,%f]", _nTotF*halfFracSigF, _nTotF*maxFracSigF));
                _work->factory(TString::Format("nBkgF[%f,%f,%f]", _nTotF*(minFracBkgF+halfFracSigF), _nTotF*minFracBkgF, _nTotF*1.5));          
            } else {
                _work->factory(TString::Format("nSigF[%f,0.5,%f]", _nTotF*0.9, _nTotF*1.5));
                _work->factory(TString::Format("nBkgF[%f,0.5,%f]", _nTotF*0.1, _nTotF*1.5));          
            }
        }
        RooFFTConvPdf* convFail = (RooFFTConvPdf*) _work->factory("FCONV::sigFail(x, sigPhysFail, sigResFail)");
        convFail->setBufferFraction(0.5);
        _work->factory("SUM::pdfFail(nSigF*sigFail, nBkgF*bkgFail)");
        
    }

    if (_work->pdf("bkgFailBackup") != nullptr) _work->factory("SUM::pdfFailBackup(nSigF*sigFail, nBkgF*bkgFailBackup)");

    if (_isMC && _work->pdf("bkgFailMC") != nullptr) {
        _work->factory("SUM::pdfFailMC(nSigF*sigFail, nBkgF*bkgFailMC)");
        _hasShape_bkgFailMC = true;
    }

}

RooFitResult* tnpFitter::manageFit(bool isPass, int attempt = 0, std::string* lastNamePDF = nullptr, double* chi2value = nullptr) {
    
    std::string pdfName = "pdfPass";
    std::string hName = "hPass";
    std::string constrainName = "constrainP";
    std::string sigPar = "nSigP";
    std::string bkgPar = "nBkgP";

    if (!isPass) {
        if (_isMC && _hasShape_bkgFailMC) {
            pdfName = (attempt==2) ? "pdfFailBackup" : "pdfFailMC";
        } else {
            pdfName = (attempt==2) ? "pdfFailBackup" : "pdfFail";
        }
        hName = "hFail";
        constrainName = "constrainF";
        sigPar = "nSigF";
        bkgPar = "nBkgF";
    }

    if (lastNamePDF != nullptr) *lastNamePDF = pdfName;

    // should check the parameters of the actual pdf being used rather than assuming that there are no constraints when attempt == 2
    //const RooArgSet* constraint = (attempt == 2) ? nullptr : _work->set(constrainName.c_str());
    const RooArgSet* constraint = _work->set(constrainName.c_str());

    if (_isMC && _hasShape_bkgFailMC) constraint = nullptr;

    RooAbsPdf *pdf = _work->pdf(pdfName.c_str());
    RooAbsData* dh =  _work->data(hName.c_str());
    RooFitResult* res = pdf->fitTo(*dh,
                                   Range("fitMassRange"),
                                   Minimizer("Minuit2"),
                                   EvalBackend("legacy"),
                                   Strategy(isPass ? _strategyPassFit : _strategyFailFit),
                                   (constraint != nullptr) ? ExternalConstraints(*constraint) : RooCmdArg::none(),
                                   SumW2Error(kFALSE), // default is false, but needs it explicitly for MC otherwise roofit complains
                                   Minos(kFALSE),
                                   Offset("bin"),
                                   PrintLevel(_printLevel),
                                   Save()
                                   );
    
    /*
    RooAbsReal * chi2 = pdf->createChi2(*((RooDataHist*) dh), Range(_xFitMin,_xFitMax));
    *chi2value = chi2->getVal();
    */

    std::string parHist = isPass ? "paramHistP" : "paramHistF";
    int auxMeas = (_work->function(parHist.c_str()) != nullptr) ? _nFitBins : 0;
    int nFloatPars = res->floatParsFinal().getSize() - auxMeas;
    
    evalChi2(pdfName, hName, nFloatPars, isPass);

    *chi2value = (isPass) ? _chi2P : _chi2F;
    int& ndof = (isPass) ? _ndofP : _ndofF;
    
    double chi2sigma = std::sqrt(2.0*ndof);
    bool goodChi2 = std::fabs(*chi2value-(double)ndof) < (100.0*chi2sigma); 

    if (attempt > 0) return res;

    if (_isMC) {
        if (goodChi2 && (res->status() == 0 or res->status() == 1)) return res;
    } 
    else {
        if (goodChi2 && res->covQual() == 3 && (res->status() == 0 or res->status() == 1)) return res;
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
    } 
    else {
        if (isPass) return res;
        else        return manageFit(isPass, 2, lastNamePDF, chi2value);
    }
    
}


int tnpFitter::fits(const std::string& title) {

    /// FC: seems to be better to change the actual range than using a fitRange in the fit itself (???)
    /// FC: I don't know why but the integral is done over the full range in the fit not on the reduced range
    _work->var("x")->setRange(_xFitMin, _xFitMax);
    _work->var("x")->setRange("fitMassRange", _xFitMin, _xFitMax);

    // TODO: check that all parameters exists
    for (auto i=_constraints.begin(); i!=_constraints.end(); i++) {
        _work->defineSet((i->first).c_str(), (i->second).c_str());
    }

    std::string lastNamePassPDF = "";
    double chi2valuePass = 0.0;
    RooFitResult* resPass = manageFit(true, 0, &lastNamePassPDF, &chi2valuePass);

    std::string lastNameFailPDF = "";
    double chi2valueFail = 0.0;
    RooFitResult* resFail = manageFit(false, 0, &lastNameFailPDF, &chi2valueFail);

    std::string bkgNamePass = "bkgPass";
    std::string bkgNameFail = "bkgFail";
    if (lastNamePassPDF.find("Backup") != std::string::npos) bkgNamePass += "Backup";
    if (lastNameFailPDF.find("Backup") != std::string::npos) {
        bkgNameFail += "Backup";
    } else {
        if (_isMC && _hasShape_bkgFailMC) bkgNameFail = "bkgFailMC"; 
    }
        
    RooPlot *pPass = _work->var("x")->frame(_xFitMin, _xFitMax);
    RooPlot *pFail = _work->var("x")->frame(_xFitMin, _xFitMax);
    pPass->SetTitle("passing probe");
    pFail->SetTitle("failing probe");

    
    _work->data("hPass")->plotOn(pPass, Name("data_pass"), MarkerSize(0.5), MarkerStyle(20));
    _work->pdf(lastNamePassPDF.c_str())->plotOn(pPass, Name("model_pass"), LineColor(kRed),  LineWidth(3));
    _work->pdf(lastNamePassPDF.c_str())->plotOn(pPass, Name("bkg_pass"),   LineColor(kBlue), LineWidth(3), Components(bkgNamePass.c_str()), LineStyle(7));
    _work->data("hPass")->plotOn(pPass, Name("data_pass"), MarkerSize(0.5), MarkerStyle(20));

    _work->data("hFail")->plotOn(pFail, Name("data_fail"), MarkerSize(0.5), MarkerStyle(20));
    _work->pdf(lastNameFailPDF.c_str())->plotOn(pFail, Name("model_fail"), LineColor(kRed),  LineWidth(3));
    _work->pdf(lastNameFailPDF.c_str())->plotOn(pFail, Name("bkg_fail"),   LineColor(kBlue), LineWidth(3), Components(bkgNameFail.c_str()), LineStyle(7));
    _work->data("hFail")->plotOn(pFail, Name("data_fail"), MarkerSize(0.5), MarkerStyle(20));

    std::string canvasName = _histname_base + "_Canv"; // TString::Format("%s_Canv",_histname_base.c_str()); 
    TCanvas * c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1150, 500);
    c->Divide(3,1);
    TPad *padText = (TPad*)c->GetPad(1);
    textParForCanvas(padText, resPass, resFail, chi2valuePass, chi2valueFail);
    c->cd(2);
    pPass->Draw();
    c->cd(3);
    pFail->Draw();

    c->cd(0);
    c->SaveAs(TString::Format("%s%s.png", _outPlotPath.c_str(), canvasName.c_str()));

    _fOut->cd();
    c->Write(canvasName.c_str(), TObject::kOverwrite);
    resPass->Write(TString::Format("%s_resP", _histname_base.c_str()), TObject::kOverwrite);
    resFail->Write(TString::Format("%s_resF", _histname_base.c_str()), TObject::kOverwrite);

    return 1;
}

double tnpFitter::getEfficiencyUncertainty(double nP, double nF, double e_nP, double e_nF) {
    double nTot = nP + nF; 
    return 1./(nTot*nTot) * std::sqrt( nP*nP* e_nF*e_nF + nF*nF * e_nP*e_nP );
}

/////// Stupid parameter dumper /////////
void tnpFitter::textParForCanvas(TPad *p, RooFitResult *resP, RooFitResult *resF, double& chi2valuePass, double& chi2valueFail) {

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

    double e_eff_corr = e_eff;
    if (!_isMC) {
        double e_nP_corr = std::max(e_nP, std::sqrt(nP));
        double e_nF_corr = std::max(e_nF, std::sqrt(nF));
        // std::cout << "Corrected stat uncertainties on nP and nF --> " << e_nP_corr << ", " << e_nF_corr << std::endl;
        e_eff_corr = getEfficiencyUncertainty(nP, nF, e_nP_corr, e_nF_corr); 
    }

    TPaveText *text1 = new TPaveText(0, 0.76, 1, 1);
    text1->SetFillColor(0);
    text1->SetBorderSize(0);
    text1->SetTextAlign(12);

    text1->AddText(TString::Format("Fit status:  pass %d, fail %d", resP->status(),  resF->status()));
    text1->AddText(TString::Format("Cov quality: pass %d, fail %d", resP->covQual(), resF->covQual()));
    //int ndofP = _nFitBins - resP->floatParsFinal().getSize() + ((_work->function("paramHistP") != nullptr) ? _nFitBins : 0);
    //int ndofF = _nFitBins - resF->floatParsFinal().getSize() + ((_work->function("paramHistF") != nullptr) ? _nFitBins : 0);
    double chi2probPass = 100.0 * TMath::Prob(chi2valuePass, _ndofP);
    double chi2probFail = 100.0 * TMath::Prob(chi2valueFail, _ndofF);
    text1->AddText(TString::Format("#Chi^{2} (prob): P %.1f/%d (%.1f%%), F %.1f/%d (%.1f%%)", chi2valuePass, _ndofP, chi2probPass, chi2valueFail, _ndofF, chi2probFail));
    //text1->SetTextFont(62);
    if (!_isMC && (e_eff_corr > e_eff) ) {
        text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f (%1.4f)",eff, e_eff, e_eff_corr));
    } else {
        text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f", eff, e_eff));
    }

    TPaveText *text = new TPaveText(0, 0, 1, 0.76);
    text->SetFillColor(0);
    text->SetBorderSize(0);
    text->SetTextAlign(12);
    std::string bkgWarning = "";
    if (_work->var("nBkgP")->getVal() <= 0.0) bkgWarning += "  nBkgP=0";
    if (_work->var("nBkgF")->getVal() <= 0.0) bkgWarning += "  nBkgF=0";

    text->AddText(TString::Format("    --- parameters %s", bkgWarning.c_str()) );

    RooArgList listParFinalP = resP->floatParsFinal();
    for (int ip=0; ip<listParFinalP.getSize(); ip++) {
        TString vName = listParFinalP[ip].GetName();
        if (!vName.Contains("_gamma_")) {
            text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f", vName.Data(), _work->var(vName)->getVal(), _work->var(vName)->getError()));
        }
    }

    RooArgList listParFinalF = resF->floatParsFinal();
    for(int ip=0; ip<listParFinalF.getSize(); ip++) {
        TString vName = listParFinalF[ip].GetName();
        if (!vName.Contains("_gamma_")) {
            text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f", vName.Data(), _work->var(vName)->getVal(), _work->var(vName)->getError()));
        }
    }

    p->cd();
    text1->Draw();
    text->Draw();

}


#endif
