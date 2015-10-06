import numpy as np
import matplotlib.pyplot as plt
import timeit as tm

n = 9
center = (n - 1) / 2

from pymks.tools import draw_microstructures
from pymks.datasets import make_delta_microstructures

X_delta = make_delta_microstructures(n_phases=2, size=(n, n, n))
draw_microstructures(X_delta[:, center])

from pymks.datasets import make_elastic_FE_strain_delta
from pymks.tools import draw_microstructure_strain

elastic_modulus = (80, 120)
poissons_ratio = (0.3, 0.3)
macro_strain = 0.02
size = (n, n, n)

t = tm.time.time()
X_delta, strains_delta = make_elastic_FE_strain_delta(elastic_modulus=elastic_modulus,
                                                      poissons_ratio=poissons_ratio,
                                                      size=size, macro_strain=macro_strain)
print 'Elapsed Time',tm.time.time() - t, 'Seconds'

draw_microstructure_strain(X_delta[0, center, :, :], strains_delta[0, center, :, :])
from pymks import MKSLocalizationModel
from pymks.bases import PrimitiveBasis

prim_basis = PrimitiveBasis(n_states=2)
model = MKSLocalizationModel(basis=prim_basis)

model.fit(X_delta, strains_delta)

from pymks.tools import draw_coeff

coeff = model.coeff
draw_coeff(coeff[center])

from pymks.datasets import make_elastic_FE_strain_random

np.random.seed(99)
t = tm.time.time()
X, strain = make_elastic_FE_strain_random(n_samples=1, elastic_modulus=elastic_modulus,
                                          poissons_ratio=poissons_ratio, size=size, macro_strain=macro_strain)
                 
print 'Elapsed Time',(tm.time.time() - t), 'Seconds'
draw_microstructure_strain(X[0, center] , strain[0, center])

t = tm.time.time()
print("The type of data that must be input to the predict function is of type: %s" % str(type(X)))
print("shape: %s" % str(X.shape))
strain_pred = model.predict(X)
print 'Elapsed Time',tm.time.time() - t,'Seconds'

from pymks.tools import draw_strains_compare

draw_strains_compare(strain[0, center], strain_pred[0, center])

from pymks.tools import draw_differences

draw_differences([strain[0, center] - strain_pred[0, center]], ['Finite Element - MKS'])