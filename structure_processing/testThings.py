import ReadMSFunction as RMS
import MKStemp as mks
import PredictPlasticity as PP

ms = RMS.readDirectory("C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/testing_data/0")
ms = ms[-1]
strains = mks.TrainPredict(False, ["C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/training_data/", "C:/Users/pkern3/Documents/MIC-AL7075-PARTICLES/testing_data/0"])
strains = strains[:,-1,...]
fips = PP.predFIPs(ms, strains)