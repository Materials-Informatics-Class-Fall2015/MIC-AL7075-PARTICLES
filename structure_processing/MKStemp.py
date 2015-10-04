import numpy as np
import os
import matplotlib.pyplot as plt
from pymks.tools import draw_microstructures
from pymks.datasets import make_delta_microstructures
from pymks.datasets import make_elastic_FE_strain_delta
from pymks.tools import draw_microstructure_strain
from pymks import MKSLocalizationModel
from pymks.bases import PrimitiveBasis
from pymks.tools import draw_coeff
from pymks.datasets import make_elastic_FE_strain_random
from pymks.tools import draw_strains_compare
from pymks.tools import draw_differences
import ReadResponses as RR
import ReadMSFunction as RMS

if __name__ == "__main__":
    import sys
    dir_test = os.getcwd()
    if(os.path.isdir(sys.argv[-1])):
        dir_test = sys.argv[-1]
  	dir_train = os.getcwd()
    if(os.path.isdir(sys.argv[-2])):
        dir_train = sys.argv[-2]	

	n = 21
	center = (n - 1) / 2
	ms_delta = make_delta_microstructures(n_phases = 2, size = (n,n,n))
	draw_microstructures(ms_delta[:,center])
	prim_basis = PrimitiveBasis(n_states = 2)
	model = MKSLocalizationModel(basis = prim_basis)
	strain_list_train = RR.readDirectory(dir_train)
	model.fit(ms_delta,strain_list_train)
	coeff = model.coeff
	draw_coeff(coeff[center])
	strain_pred = model.predict(ms_list)
	draw_strains_compare(strain_list_test[0,center],strain_pred[0,center])
	draw_differences([strain_list_test[0,center]-strain_pred[0,center]],['Strain from Testing - MKS'])
