import math
from array import array
import ROOT

class efficiency:

    iAltBkgModel = 0
    iAltSigModel = 1
    iAltMCSignal = 2
    iAltTagSelec = 3
    iPUup        = 4
    iPUdown      = 5
    iAltFitRange = 6

    def __init__(self, abin):
        self.ptBin   = abin
        self.effData = -1
        self.effMC   = -1
        self.altEff  = [-1]*7
        self.syst    = [-1]*7


    def __init__(self, ptBin, etaBin, effData, errEffData, effMC, errEffMC, 
                 effAltSigModel=-1, errAltSigModel=-1, effAltBkgModel=-1, errAltBkgModel=-1,
                 effMCAltSigModel=-1, errMCAltSigModel=-1, effMCAltBkgModel=-1, errMCAltBkgModel=-1,
                 effAltTagSel=-1):
        """
        """
        self.ptBin      = ptBin
        self.etaBin     = etaBin
        self.effData    = effData
        self.effMC      = effMC
        self.errEffData = errEffData
        self.errEffMC   = errEffMC
        self.effAltSig = effAltSigModel
        self.errAltSig = errAltSigModel
        self.effAltBkg = effAltBkgModel
        self.errAltBkg = errAltBkgModel
        self.effMCAltSig = effMCAltSigModel
        self.errMCAltSig = errMCAltSigModel
        self.effMCAltBkg = effMCAltBkgModel
        self.errMCAltBkg = errMCAltBkgModel
        self.altEff = {
            "altSig": [effAltSigModel, errAltSigModel],
            "altBkg": [effAltBkgModel, errAltBkgModel],
            "MC_altSig": [effMCAltSigModel, errMCAltSigModel],
            "MC_altBkg": [effMCAltBkgModel, errMCAltBkgModel],
            "tagSel": [effAltTagSel, -1]
        }
        self.syst   = {}


    def __str__(self):
        """
        """
        strout = f'{self.etaBin[0]:2.3f}\t{self.etaBin[1]:2.3f}\t{self.ptBin[0]:2.1f}\t{self.ptBin[1]:2.1f}'
        strout += f'\t{self.effData:2.4f}\t{self.errEffData:2.4f}\t{self.effMC:2.4f}\t{self.errEffMC:2.4f}'
        
        for k, v in self.altEff.items():
            if v[0] > 0:
                strout += f'\t{v[0]:2.4f}\t{v[1]:2.4f}'
        return strout


    def __add__(self, eff):
        """
        """

        ## TO BE CHECKED

        if self.effData < 0 :
            return eff.deepcopy()
        if eff.effData < 0 :
            return self.deepcopy()
        
        ptbin  = self.ptBin
        etabin = self.etaBin
        errEffData = self.errEffData if self.errEffData else 1.
        efferrEffData = eff.errEffData if eff.errEffData else 1.
        
        errData2 = 1.0 / (1.0/(errEffData*errEffData)+1.0/(efferrEffData*efferrEffData))
        wData1   = 1.0 / (errEffData * errEffData) * errData2
        wData2   = 1.0 / (efferrEffData * efferrEffData) * errData2
        newEffData      = wData1 * self.effData + wData2 * eff.effData
        newErrEffData   = math.sqrt(errData2)
        
        #        errMC2 = 1.0 / (1.0/(self.errEffMC*self.errEffMC)+1.0/(eff.errEffMC*eff.errEffMC))
        #wMC1   = 1.0 / (self.errEffMC * self.errEffMC) * errMC2
        #wMC2   = 1.0 / (eff .errEffMC * eff .errEffMC) * errMC2
        newEffMC      = wData1 * self.effMC + wData2 * eff.effMC
        newErrEffMC   = 0.00001#math.sqrt(errMC2)

        newEffAltBkgModel = wData1 * self.altEff[self.iAltBkgModel] + wData2 * eff.altEff[self.iAltBkgModel]
        newEffAltSigModel = wData1 * self.altEff[self.iAltSigModel] + wData2 * eff.altEff[self.iAltSigModel]
        newEffAltMCSignal = wData1 * self.altEff[self.iAltMCSignal] + wData2 * eff.altEff[self.iAltMCSignal]
        newEffAltTagSelec = wData1 * self.altEff[self.iAltTagSelec] + wData2 * eff.altEff[self.iAltTagSelec]

        #effout = efficiency(ptbin, etabin, newEffData, newErrEffData, newEffMC, newErrEffMC, newEffAltBkgModel,newEffAltSigModel,newEffAltMCSignal,newEffAltTagSelec)
        return None
    

    @staticmethod
    def getSystematicNames():
        return ['statData', 'statMC', 'altSignalModel', 'altBkgModel', 'altMCEff', 'altTagSelection']


    def combineSyst(self):
        """
        """
        systAltSig = max(self.altEff["altSig"][0], 0)
        systAltBkg = max(self.altEff["altBkg"][0], 0)
        systMC_AltSig = max(self.altEff["MC_altSig"][0], 0)
        systMC_AltBkg = max(self.altEff["MC_altBkg"][0], 0)

        self.syst[0] = self.errEffData
        self.syst[1] = systAltBkg
        self.syst[2] = systAltSig
        self.syst[3] = systMC_AltSig
        self.syst[4] = systMC_AltBkg

        for i in range(len(self.syst)):
            self.syst[i] = self.syst[i]/self.effData if self.effData > 0 else 1.
            if i==0:
                self.syst[i] = 1. # stat error, not included in this systematic combination
         
        self.systCombinedVar = 0  # combined systematic variation
        for s in self.syst:
            s_add_upVar = max(s, 2-s)
            self.systCombinedVar += s_add_upVar * s_add_upVar

        self.systCombinedVar = (2-math.sqrt(self.systCombinedVar), math.sqrt(self.systCombinedVar))  # variation down-up

    
class efficiencyManager: 



    def __init__(self):
        self.effList = {}

    def __str__(self):
        outStr = ''
        for ptBin in self.effList.keys():
            for etaBin in self.effList[ptBin].keys():
                outStr += f"{str(self.effList[ptBin][etaBin])}\n"
        return outStr

    def addEfficiency(self, eff):
        """
        """
        if not eff.ptBin in self.effList:
            self.effList[eff.ptBin] = {}
        self.effList[eff.ptBin][eff.etaBin] = eff

    def combineSyst(self):
        """
        """
        for ptBin in self.effList.keys():
            for etaBin in self.effList[ptBin].keys():
                self.effList[ptBin][etaBin].combineSyst()
                self.effList[ptBin][etaBin].combineSyst()
                                   
    def symmetrizeSystVsEta(self):
        """
        """
        for ptBin in self.effList.keys():
            for etaBin in self.effList[ptBin].keys():
                if etaBin[0] >= 0 and etaBin[1] > 0:
                    etaBinPlus  = etaBin
                    etaBinMinus = (-etaBin[1],-etaBin[0])
                    
                    effPlus  = self.effList[ptBin][etaBinPlus]
                    effMinus = None
                    if etaBinMinus in self.effList[ptBin]: #python3self.effList[ptBin].has_key(etaBinMinus):
                        effMinus =  self.effList[ptBin][etaBinMinus] 

                    if effMinus is None:
                        self.effList[ptBin][etaBinMinus] = effPlus
                        pass # print(" ---- efficiencyList: I did not find -eta bin!!!")
                    else:
                        #### fix statistical errors if needed
                        if    effPlus.errEffData <= 0.00001 and effMinus.errEffData > 0.00001: 
                            self.effList[ptBin][etaBinPlus ].errEffData = effMinus.errEffData
                        elif effMinus.errEffData <= 0.00001 and effPlus .errEffData > 0.00001: 
                            self.effList[ptBin][etaBinMinus].errEffData = effPlus.errEffData
                        else:
                            self.effList[ptBin][etaBinPlus ].errEffData = (effMinus.errEffData+effPlus.errEffData)/2.
                            self.effList[ptBin][etaBinMinus].errEffData = (effMinus.errEffData+effPlus.errEffData)/2.

                        if   effPlus.errEffMC <= 0.00001 and effMinus.errEffMC > 0.00001: 
                            self.effList[ptBin][etaBinPlus ].errEffMC = effMinus.errEffMC
                        elif effMinus.errEffMC <= 0.00001 and effPlus.errEffMC > 0.00001: 
                            self.effList[ptBin][etaBinMinus].errEffMC = effPlus.errEffMC
                        else:
                            self.effList[ptBin][etaBinPlus ].errEffMC = (effMinus.errEffMC+effPlus.errEffMC)/2.
                            self.effList[ptBin][etaBinMinus].errEffMC = (effMinus.errEffMC+effPlus.errEffMC)/2.

                            
                        for isyst in range(4):
                            if abs(effPlus.altEff[isyst] - effMinus.altEff[isyst]) < 0.10:
                                averageSyst = (effPlus.altEff[isyst] +  effMinus.altEff[isyst]) / 2
                                self.effList[ptBin][etaBinPlus ].altEff[isyst] = averageSyst
                                self.effList[ptBin][etaBinMinus].altEff[isyst] = averageSyst
                            else:
                                averageSyst = (effPlus.altEff[isyst] +  effMinus.altEff[isyst]) / 2
                                print("issue, I am averaging but the efficiencies are quite different in 2 etaBins")
                                print(" --- syst: ", isyst)
                                print(str(self.effList[ptBin][etaBinPlus ]))
                                print(str(self.effList[ptBin][etaBinMinus]))
                                print("   eff[+] = ",  self.effList[ptBin][etaBinPlus ].altEff[isyst])
                                print("   eff[-] = ",  self.effList[ptBin][etaBinMinus].altEff[isyst])
                                self.effList[ptBin][etaBinPlus ].altEff[isyst] = averageSyst
                                self.effList[ptBin][etaBinMinus].altEff[isyst] = averageSyst

    def get1DGraphList(self, var, doScaleFactor=False, effMC=False, typeErr=""):
        """
        """
        listOfGraphs = {}

        if var not in ["eta", "pt"]:
            print(f" --- efficiencyManager: var {var} not found")
            return listOfGraphs
        
        val_attr = "effData" if not effMC else "effMC"
        err_attr = "errEffData" if not effMC else "errEffMC"

        

        if var == "eta":
            listOfGraphs = {ptBin: [] for ptBin in self.effList.keys()}
        else:
            listOfGraphs = {etaBin: [] for ptBin in self.effList.keys() for etaBin in self.effList[ptBin].keys()}

        for ptBin in self.effList.keys():

            for etaBin in self.effList[ptBin].keys():

                val = getattr(self.effList[ptBin][etaBin], val_attr)

                if typeErr in ["", "stat"]:
                    err_val = getattr(self.effList[ptBin][etaBin], err_attr)
                elif typeErr == "allSyst":
                    self.combineSyst()
                    _, var_up = self.effList[ptBin][etaBin].systCombined 
                    err_val = val * (var_up-1)
                else:
                    if not typeErr in self.altEff:
                        print(f" --- efficiencyManager: {typeErr} not found in altEff")
                        continue
                    err_val = self.effList[ptBin][etaBin].altEff[typeErr][0] / val
                
                if doScaleFactor:
                    val /= self.effList[ptBin][etaBin].effMC
                    err_val /= self.effList[ptBin][etaBin].effMC  # uncertainty of MC efficiency is supposed to be small

                if var == "eta":
                    listOfGraphs[ptBin].append({
                        'min': etaBin[0], 'max': etaBin[1], 'val': val, 'err': err_val
                    })
                else:
                    listOfGraphs[etaBin].append({
                        'min': ptBin[0], 'max': ptBin[1], 'val': val, 'err': err_val
                    })
                
                var_differential, var_analysis = (ptBin, etaBin) if var=="eta" else (etaBin, ptBin)

                listOfGraphs[var_differential].append({
                    'min': var_analysis[0], 'max': var_analysis[1], 'val': val, 'err': err_val
                    })

                                                  
        return listOfGraphs
    
    
    def ptEtaScaleFactor_2DHisto(self, typePlot="eff", isMC=False, typeEff="nominal"):
        """
        """

        xbins = sorted({edge for ptBin in self.effList for etaBin in self.effList[ptBin] for edge in etaBin})
        ybins = sorted({edge for ptBin in self.effList for edge in ptBin})
        xbinsTab = array('d', xbins)
        ybinsTab = array('d', ybins)

        plot_types = {
            "eff": ("h2_eff", "Lepton efficiencies "),
            "sf": ("h2_sf", "Lepton scale factors "),
            "err": ("h2_statUnc", "Lepton eff. stat. unc. "),
            "syst": ("h2_systUnc", "Lepton eff. syst. unc. ")
        }
        isMC_types = {
            False: ("Data", "data "),
            True: ("MC", "MC ")
        }
        eff_types = {
            "nominal": ("Nominal", "Nominal"),
            "altSig": ("AltSig", "AltSig"),
            "altBkg": ("AltBkg", "AltBkg"),
            "allSyst": ("", "")
        }
        baseName, baseTitle = plot_types[typePlot]
        effName, effTitle = eff_types[typeEff]

        if typePlot not in plot_types.keys():
            print(" --- efficiencyManager: typePlot not found")
            return None

        if typeEff not in eff_types.keys():
            print(" --- efficiencyManager: typeEff not found")
            return None
        
        if (typePlot=="syst" and typeEff=="nominal") or (typePlot=="sf" and isMC):
            print(f" --- efficiencyManager: typePlot={typePlot} not compatible with selected settings")
            return None 

        hname = baseName + isMC_types[isMC][0] + effName
        htitle = baseTitle + isMC_types[isMC][1] + effTitle
        
        # names of the attributes of the efficiency class to be used as value and error
        var_attribute = hname.replace("h2_", "").replace("sf", "eff").replace("Data", "" if typeEff!="nominal" else "Data").replace("Nominal", "")
        err_attribute = var_attribute.replace("eff", "err" if typeEff!="nominal" else "errEff")
        
        if typePlot == "err":
            var_attribute = "err"
            var_attribute += "Eff" if typeEff == "nominal" else ""
            var_attribute += "MC" if isMC else ""
            var_attribute += effName
            var_attribute = var_attribute.replace("Nominal", "Data" if not isMC else "")
            err_attribute = None
        elif typePlot == "syst":
            if typeEff == "allSyst":
                self.combineSyst()
                var_attribute = "systCombinedVar"
            else:
                var_attribute = var_attribute.replace("systUnc", "eff")
            err_attribute = None

        print(hname, htitle, "   ", var_attribute, err_attribute)

        h2 = ROOT.TH2D(hname, htitle, len(xbinsTab)-1, xbinsTab, len(ybinsTab)-1, ybinsTab)
        h2.Sumw2()

        ## init histogram efficiencies and errors to 100%
        for ix in range(1, h2.GetXaxis().GetNbins()+1):
            for iy in range(1, h2.GetYaxis().GetNbins()+1):
                h2.SetBinContent(ix, iy, 1)
                h2.SetBinError  (ix, iy, 1)

        for ix in range(1, h2.GetXaxis().GetNbins()+1):
            for iy in range(1, h2.GetYaxis().GetNbins()+1):

                for ptBin in self.effList.keys():
                    if h2.GetYaxis().GetBinLowEdge(iy) < ptBin[0] or h2.GetYaxis().GetBinUpEdge(iy) > ptBin[1]:
                        continue
                    for etaBin in self.effList[ptBin].keys():
                        if h2.GetXaxis().GetBinLowEdge(ix) < etaBin[0] or h2.GetXaxis().GetBinUpEdge(ix) > etaBin[1]:
                            continue
                        
                        #print(var_attribute, err_attribute)
                        val = getattr(self.effList[ptBin][etaBin], var_attribute)
                        err = getattr(self.effList[ptBin][etaBin], err_attribute) if err_attribute else 0

                        if typePlot == "sf":
                            val /= self.effList[ptBin][etaBin].effMC
                            err /= self.effList[ptBin][etaBin].effMC
                        elif typePlot == "syst":
                            if typeEff == "allSyst":
                                val = val[0]
                            else:
                                val /= self.effList[ptBin][etaBin].effData
                                err /= self.effList[ptBin][etaBin].effData

                        h2.SetBinContent(ix, iy, val)
                        h2.SetBinError  (ix, iy, err)

        h2.GetXaxis().SetTitle("#eta")
        h2.GetYaxis().SetTitle("p_{T} (GeV)")
        return h2
        
                                
    



def makeTGraphFromListEff(listOfEfficiencies, keyMin, keyMax):
    """
    """
    grOut = ROOT.TGraphErrors(len(listOfEfficiencies))
    
    for ip, point in enumerate(listOfEfficiencies):
        grOut.SetPoint(     ip, (point[keyMin]+point[keyMax])/2., point['val'])
        grOut.SetPointError(ip, (point[keyMax]-point[keyMin])/2., point['err'])

    return grOut
    
