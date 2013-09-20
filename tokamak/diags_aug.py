class DCN:
    def __init__(self, shn, eqi=None):
        from digitizer_aug import DigitizerAUGDCR
        from diaggeom_aug import DCNGeom

        self.digitizer = DigitizerAUGDCR(shn)
        self.geom = DCNGeom()
        self.eqi = eqi
        

