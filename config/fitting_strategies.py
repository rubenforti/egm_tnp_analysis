from libPython import fitUtils_reco_trk as fitUtils

fitStrategies = {
    "reco" : {
        "nominal" : fitUtils.histFitterNominal,
        "altSig"  : fitUtils.histFitterAltSig,
        "altBkg"  : fitUtils.histFitterAltBkg
    },
    "tracking" : {
        "nominal" : fitUtils.histFitterNominal,
        "altSig"  : fitUtils.histFitterAltSig,
        "altBkg"  : fitUtils.histFitterAltBkg
    },
    "idip" : {
        "nominal" : fitUtils.histFitterAltBkg,
        "altSig"  : fitUtils.histFitterAllAnalytic,
        "altBkg"  : None
    },
    "trigger" : {
        "nominal" : fitUtils.histFitterAltBkg,
        "altSig"  : fitUtils.histFitterAllAnalytic,
        "altBkg"  : None
    },
    "iso" : {
        "nominal" : fitUtils.histFitterAltBkg,
        "altSig"  : fitUtils.histFitterAllAnalytic,
        "altBkg"  : None
    },
    "veto" : {
        "nominal" : fitUtils.histFitterAltBkg,
        "altSig"  : fitUtils.histFitterAllAnalytic,
        "altBkg"  : None
    }
}