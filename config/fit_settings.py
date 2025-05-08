"""
"""

def fitBinning(typeflag, extendLowPtForVeto=False):
    """
    """

    binning_eta = [round(-2.4+0.1*i, 2) for i in range(49)]

    if typeflag == 'reco':
        massbins, massmin, massmax = 60, 60, 120
        binning_pt = [24., 26., 30., 34., 38., 42., 46., 50., 55., 60., 65.]
        lowPtBinsVeto = [10., 15., 20.]
    
    elif typeflag == 'tracking':
        massbins, massmin, massmax = 80, 50, 130
        binning_pt = [24., 35., 45., 55., 65.]
        lowPtBinsVeto = [10., 15.]
    
    else:
        massbins, massmin, massmax = 60, 60, 120
        binning_pt  = [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.]
        lowPtBinsVeto = [10., 15., 20.]

    if extendLowPtForVeto: 
        binning_pt = lowPtBinsVeto + binning_pt

    binningDef_mass = {"nbins": massbins, "min": massmin, "max": massmax}
    binningDef = {
        'eta' : {'var': 'eta', 'type': 'float', 'bins': binning_eta},
        'pt'  : {'var': 'pt',  'type': 'float', 'bins': binning_pt}
    }

    return binningDef_mass, binningDef


def fitParsAndShapes(typeflag):
    """
    """

    parsAndShapes = {
        "tnpParNominal": [],
        "tnpParAltSig": [],
        "tnpParAltBkg": [],
        
        "tnpShapesNominal": [],
        "tnpShapesAltSig": [],
        "tnpShapesAltBkg": [],

        "parConstraints": [],
    }

    ###########################################################################

    if typeflag == 'reco':

        tnpShapesNominal = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpParNominal = [
            "meanP[-0., -5., 5.]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0., -3., 3.]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
        ]

        tnpShapesAltSig = [
            "BreitWigner::sigPhysPass(x, 91.1876, 2.4952)",
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "RooCBExGaussShape::sigResPass(x, meanP, sigmaP, alphaP, nP, sigmaP_2, tailLeft)",
            "Exponential::bkgPass(x, expalphaP)",
            "BreitWigner::sigPhysFail(x, 91.1876, 2.4952)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "RooCBExGaussShape::sigResFail(x, meanF, sigmaF, alphaF, nF, sigmaF_2, tailLeft)",
        ]
        tnpParAltSig = [
            "tailLeft[-1]", "tailLeft[1]",
            "meanP[-0., -5., 5.]", "sigmaP[1.0, 0.7, 6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3.0, 0.01, 5.0]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[-0., -5., 5.]", "sigmaF[2.0, 0.7, 5.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3.0, 0.1,  5.0]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]",
        ]

        tnpShapesAltBkg = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
            "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F})",
        ]
        tnpParAltBkg = [
            "meanP[-0., -5., 5.]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0., -3., 3.]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
            "acmsF[60., 40., 130.]", "betaF[5.0, 0.1, 40.0]", "gammaF[0.1, 0.0, 1.0]", "peakF[90.0]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]"
        ]
        
        
        parConstraints = [
            # Passing
            #"Gaussian::constrainP_acmsP(acmsP,90,50)",
            #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
            #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
            # Failing
            "Gaussian::constrainF_acmsF(acmsF, 90, 50)",
            "Gaussian::constrainF_betaF(betaF, 5.0, 25.0)",
            "Gaussian::constrainF_gammaF(gammaF, 0.5, 0.8)",
        ]

    ###########################################################################

    elif typeflag == 'tracking':

        tnpShapesNominal = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpParNominal = [
            "meanP[-0., -5., 5.]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0., -3., 3.]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
        ]

        tnpShapesAltSig = [
            "BreitWigner::sigPhysPass(x, 91.1876, 2.4952)",
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "RooCBExGaussShape::sigResPass(x, meanP, sigmaP, alphaP, nP, sigmaP_2, tailLeft)",
            "BreitWigner::sigPhysFail(x, 91.1876, 2.4952)",
            "Exponential::bkgPass(x, expalphaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "RooCBExGaussShape::sigResFail(x, meanF, sigmaF, alphaF, nF, sigmaF_2, tailLeft)",
        ]
        tnpParAltSig = [
            "tailLeft[1]", "tailLeft[-1]", 
            "meanP[-0.,  -5.,  5.]", "sigmaP[1.0, 0.7,  6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3.0, 0.01, 5.0]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[-0., -12., 12.]", "sigmaF[2.0, 0.7, 12.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3.0, 0.01, 5.0]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]"
        ]
        tnpParAltSig_trackingHighPt = [
            "tailLeft[1]", "tailLeft[-1]", 
            "meanP[0., -5.,  5.]", "sigmaP[1.0, 0.7,  6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3., 0., 5.]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[4., -1., 15.]", "sigmaF[2.0, 0.7, 15.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3., 0., 5.]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]"
        ]

        tnpShapesAltBkg = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
            "Chebychev::bkgFailBackup(x, {c1F,c2F,c3F,c4F})",
        ]
        tnpParAltBkg = [
            "meanP[-0., -5., 5.]", "sigmaP[0.5, 0.1,  5.0]",
            "meanF[-0., -5., 5.]", "sigmaF[0.5, 0.02, 3.0]",
            "expalphaP[0., -5., 5.]",
            "acmsF[60., 40., 130.]", "betaF[5.0, 0.1, 40.0]", "gammaF[0.1, 0.0, 1.0]", "peakF[90.0]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]", "c4F[-0.5, -1.0, 1.0]"
        ]

        parConstraints = [
            # Passing
            #"Gaussian::constrainP_acmsP(acmsP,90,50)",
            #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
            #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
            # Failing
            "Gaussian::constrainF_acmsF(acmsF, 90, 50)",
            "Gaussian::constrainF_betaF(betaF, 5.0, 25.0)",
            "Gaussian::constrainF_gammaF(gammaF, 0.5, 0.8)",
        ]

        # tnpParNominal.extend(["maxFracSigF[0.5]"])
        tnpParAltBkg.extend(["maxFracSigF[0.5]"])
        tnpParAltSig.extend(["maxFracSigF[0.5]"])
        tnpParAltSig_trackingHighPt.extend(["maxFracSigF[0.5]"])

    ###########################################################################

    else:

        tnpShapesNominal = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
            "Exponential::bkgFail(x, expalphaF)",
            "Chebychev::bkgFailBackup(x, {c1F,c2F,c3F})",
        ]
        tnpParNominal = [
            "meanP[-0., -5., 5.]", "sigmaP[0.5, 0.1, 5.0]",
            "meanF[-0., -5., 5.]", "sigmaF[0.5, 0.1, 5.0]",
            "expalphaP[0., -5., 5.]",
            "expalphaF[0., -5., 5.]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]",
        ]

        tnpShapesAltSig = [
            "BreitWigner::sigPhysPass(x, 91.1876, 2.4952)",
            "RooCBExGaussShape::sigResPass(x, meanP, sigmaP, alphaP, nP, sigmaP_2, tailLeft)",
            "Exponential::bkgPass(x, expalphaP)",
            "BreitWigner::sigPhysFail(x, 91.1876, 2.4952)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "RooCBExGaussShape::sigResFail(x, meanF, sigmaF, alphaF, nF, sigmaF_2, tailLeft)",
            "Gaussian::sigFsrFail(x, fsrMeanF, fsrSigmaF)",
            "Exponential::bkgFail(x, expalphaF)",
            "Chebychev::bkgFailBackup(x, {c1F,c2F,c3F})",
        ]
        tnpParAltSig = [
            "tailLeft[1]", "tailLeft[-1]", 
            "meanP[-0., -5., 5.]", "sigmaP[1.0, 0.7,  6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3.0, 0.01, 5.0]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[-0., -5., 5.]", "sigmaF[2.0, 0.7, 15.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3-0, 0.01, 5.0]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]",
            "fsrMeanF[70., 65., 80.]", "fsrSigmaF[1.0, 1.2, 5.0]",
            "expalphaF[0., -5., 5.]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]",
        ]

        tnpShapesAltBkg = [
            # Alternative background trategy not currently used for steps out of "reco" and "tracking"
        ]
        tnpParAltBkg = [
            # Alternative background trategy not currently used for steps out of "reco" and "tracking"
        ]

        parConstraints = [
            # Passing
            #"Gaussian::constrainP_acmsP(acmsP,90,50)",
            #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
            #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
            # Failing
            "Gaussian::constrainF_acmsF(acmsF, 90, 50)",
            "Gaussian::constrainF_betaF(betaF, 5.0, 25.0)",
            "Gaussian::constrainF_gammaF(gammaF, 0.5, 0.8)",
        ]
    
    ###########################################################################

    parsAndShapes["tnpParNominal"] = tnpParNominal
    parsAndShapes["tnpParAltSig"] = tnpParAltSig
    parsAndShapes["tnpParAltBkg"] = tnpParAltBkg 
    
    parsAndShapes["tnpShapesNominal"] = tnpShapesNominal
    parsAndShapes["tnpShapesAltSig"] = tnpShapesAltSig
    parsAndShapes["tnpShapesAltBkg"] = tnpShapesAltBkg
    
    parsAndShapes["parConstraints"] = parConstraints

    if typeflag == 'tracking':
        parsAndShapes["tnpParAltSig_trackingHighPt"] = tnpParAltSig_trackingHighPt


    return parsAndShapes


    