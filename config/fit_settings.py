"""
"""

def fitGeneralSettings(typeflag):
    """    
    """
    pass

def fitParsAndShapes(typeflag):
    """
    """

    parsAndShapes = {
        "tnpParNomFit": [],
        "tnpParAltSigFit": [],
        "tnpParAltBkgFit": [],
        "tnpShapesNominal": [],
        "tnpShapesAltSig": [],
        "tnpShapesAltBkg": [],
        "parConstraints": [],
    }

    if typeflag == 'reco':
        tnpParNomFit = [
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0.0, -3.0, 3.0]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
        ]
        tnpParAltSigFit = [
            "tailLeft[-1]", "tailLeft[1]",
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[1.0, 0.7, 6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3.0, 0.01, 5.0]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[-0.0, -5.0, 5.0]", "sigmaF[2.0, 0.7, 5.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3-0, 0.1,  5.0]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]"
        ]
        tnpParAltBkgFit = [
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0.0, -3.0, 3.0]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
            "acmsF[60., 40., 130.]", "betaF[5.0, 0.1, 40.0]", "gammaF[0.1, 0.0, 1.0]", "peakF[90.0]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]", "c4F[-0.5, -1.0, 1.0]"
        ]
        tnpShapesNominal = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpShapesAltSig = [
            "BreitWigner::sigPhysPass(x, 91.1876, 2.4952)",
            "RooCBExGaussShape::sigResPass(x, meanP, sigmaP, alphaP, nP, sigmaP_2, tailLeft)",
            "BreitWigner::sigPhysFail(x, 91.1876, 2.4952)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "RooCBExGaussShape::sigResFail(x, meanF, sigmaF, alphaF, nF, sigmaF_2, tailLeft)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpShapesAltBkg = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
            "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F})",
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

    elif typeflag == 'tracking':
        tnpParNomFit = [
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0.0, -3.0, 3.0]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
        ]
        tnpParAltSigFit = [
            "tailLeft[1]", "tailLeft[-1]", 
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[1.0, 0.7, 6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3.0, 0.01, 5.0]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[-0.0, -5.0, 5.0]", "sigmaF[2.0, 0.7, 5.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3-0, 0.1,  5.0]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]"
        ]
        # for pt >= 55 and tracking (se also note above)
        tnpParAltSigFitTrackingHighPt = [
            "tailLeft[1]", "tailLeft[-1]", 
            "meanP[0.0, -5.0, 5.0]",  "sigmaP[1, 0.7, 6.0]",  "alphaP[2.0, 1.2, 3.5]", "nP[3., 0., 5.]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[4.0, -1.0, 15.0]", "sigmaF[2, 0.7, 15.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3., 0., 5.]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]"
        ]
        tnpParAltBkgFit = [
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0.0, -3.0, 3.0]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
            "acmsF[60., 40., 130.]", "betaF[5.0, 0.1, 40.0]", "gammaF[0.1, 0.0, 1.0]", "peakF[90.0]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]", "c4F[-0.5, -1.0, 1.0]"
        ]
        tnpShapesNominal = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpShapesAltSig = [
            "BreitWigner::sigPhysPass(x, 91.1876, 2.4952)",
            "RooCBExGaussShape::sigResPass(x, meanP, sigmaP, alphaP, nP, sigmaP_2, tailLeft)",
            "BreitWigner::sigPhysFail(x, 91.1876, 2.4952)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "RooCBExGaussShape::sigResFail(x, meanF, sigmaF, alphaF, nF, sigmaF_2, tailLeft)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpShapesAltBkg = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
            "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F})",
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
        tnpParNomFit.extend(["maxFracSigF[0.5]"])
        tnpParAltSigFit.extend(["maxFracSigF[0.5]"])
        tnpParAltSigFitTrackingHighPt.extend(["maxFracSigF[0.5]"])
        parConstraints = [
            #"Gaussian::constrainP_acmsP(acmsP,90,50)",
            #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
            #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
            # failing
            "Gaussian::constrainF_acmsF(acmsF,90,50)",
            "Gaussian::constrainF_betaF(betaF,5.0,25.0)",
            "Gaussian::constrainF_gammaF(gammaF,0.5,0.8)",
        ]

    else:
        tnpParNomFit = [
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0.0, -3.0, 3.0]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
        ]
        tnpParAltSigFit = [
            "tailLeft[1]", "tailLeft[-1]", 
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[1.0, 0.7, 6.0]", "alphaP[2.0, 1.2, 3.5]", "nP[3.0, 0.01, 5.0]", "sigmaP_2[1.5, 0.5, 6.0]",
            "meanF[-0.0, -5.0, 5.0]", "sigmaF[2.0, 0.7, 5.0]", "alphaF[2.0, 1.2, 3.5]", "nF[3-0, 0.1,  5.0]", "sigmaF_2[2.0, 0.5, 6.0]",
            "expalphaP[0., -5., 5.]"
            "fsrMeanF[70.0, 65.0, 80.0]", "fsrSigmaF[1.0, 1.2, 5.0]"
        ]
        tnpParAltBkgFit = [
            "meanP[-0.0, -5.0, 5.0]", "sigmaP[0.5, 0.1,  3.0]",
            "meanF[-0.0, -3.0, 3.0]", "sigmaF[0.5, 0.01, 2.0]",
            "expalphaP[0., -5., 5.]",
            "acmsF[60., 40., 130.]", "betaF[5.0, 0.1, 40.0]", "gammaF[0.1, 0.0, 1.0]", "peakF[90.0]",
            "c1F[0.0, -1.0, 1.0]", "c2F[-0.5, -1.0, 1.0]", "c3F[0.0, -1.0, 1.0]", "c4F[-0.5, -1.0, 1.0]"
        ]
        tnpShapesNominal = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpShapesAltSig = [
            "BreitWigner::sigPhysPass(x, 91.1876, 2.4952)",
            "RooCBExGaussShape::sigResPass(x, meanP, sigmaP, alphaP, nP, sigmaP_2, tailLeft)",
            "BreitWigner::sigPhysFail(x, 91.1876, 2.4952)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "RooCBExGaussShape::sigResFail(x, meanF, sigmaF, alphaF, nF, sigmaF_2, tailLeft)",
            "Gaussian::sigFsrFail(x, fsrMeanF, fsrSigmaF)"
            "Exponential::bkgPass(x, expalphaP)",
        ]
        tnpShapesAltBkg = [
            "Gaussian::sigResPass(x, meanP, sigmaP)",
            "Gaussian::sigResFail(x, meanF, sigmaF)",
            "Exponential::bkgPass(x, expalphaP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
            "Chebychev::bkgFailBackup(x,{c1F,c2F,c3F})",
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
        tnpParNomFit.extend(["maxFracSigF[0.5]"])
        tnpParAltSigFit.extend(["maxFracSigF[0.5]"])
        parConstraints = [
            #"Gaussian::constrainP_acmsP(acmsP,90,50)",
            #"Gaussian::constrainP_betaP(betaP,0.05,0.25)",
            #"Gaussian::constrainP_gammaP(gammaP,0.5,0.8)",
            # failing
            "Gaussian::constrainF_acmsF(acmsF,90,50)",
            "Gaussian::constrainF_betaF(betaF,5.0,25.0)",
            "Gaussian::constrainF_gammaF(gammaF,0.5,0.8)",
        ]

    parsAndShapes["tnpParNomFit"] = tnpParNomFit
    parsAndShapes["tnpParAltSigFit"] = tnpParAltSigFit
    parsAndShapes["tnpParAltBkgFit"] = tnpParAltBkgFit
    parsAndShapes["tnpShapesNominal"] = tnpShapesNominal
    parsAndShapes["tnpShapesAltSig"] = tnpShapesAltSig
    parsAndShapes["tnpShapesAltBkg"] = tnpShapesAltBkg
    parsAndShapes["parConstraints"] = parConstraints

    if typeflag == 'tracking':
        parsAndShapes["tnpParAltSigFitTrackingHighPt"] = tnpParAltSigFitTrackingHighPt

    return parsAndShapes


    