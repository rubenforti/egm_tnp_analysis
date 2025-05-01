fitStrategies = {
    "reco" : {
        "Nominal" : "histFitterAllTemplate",
        "AltSig"  : "histFitterAnalyticSig",
        "AltBkg"  : "histFitterAnalyticBkg"
    },
    "tracking" : {
        "Nominal" : "histFitterAllTemplate",
        "AltSig"  : "histFitterAnalyticSig",
        "AltBkg"  : "histFitterAnalyticBkg"
    },
    "idip" : {
        "Nominal" : "histFitterAnalyticBkg",
        "AltSig"  : "histFitterAllAnalytic",
        "AltBkg"  : None
    },
    "trigger" : {
        "Nominal" : "histFitterAnalyticBkg",
        "AltSig"  : "histFitterAllAnalytic",
        "AltBkg"  : None
    },
    "iso" : {
        "Nominal" : "histFitterAnalyticBkg",
        "AltSig"  : "histFitterAllAnalytic",
        "AltBkg"  : None
    },
    "veto" : {
        "Nominal" : "histFitterAnalyticBkg",
        "AltSig"  : "histFitterAllAnalytic",
        "AltBkg"  : None
    }
}