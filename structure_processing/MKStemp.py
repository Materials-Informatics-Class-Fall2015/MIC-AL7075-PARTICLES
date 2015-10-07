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
        
    labels = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
    n = 21
    center = (n - 1) / 2
    
    ms_delta = make_delta_microstructures(n_phases = 2, size = (n,n,n))
    # draw_microstructures(ms_delta[:,center])
    prim_basis = PrimitiveBasis(n_states = 2)
    
    strain_list_train = RR.readDirectory(dir_train)
    tensor_comp = strain_list_train.shape[0]
    strain_list_test = RR.readDirectory(dir_test)
    print(strain_list_test.shape)
    num_ms = strain_list_test.shape[1]
    ms_list = RMS.readDirectory(dir_test)
    output_data = strain_list_test
    
    models = []
    for i in range(tensor_comp):
        temp_train = strain_list_train[i]
        model = MKSLocalizationModel(basis = prim_basis)
        model.fit(ms_delta,temp_train)
        coeff = model.coeff
        models.append(model)
        # draw_coeff(coeff[center])
    
    new_shape = strain_list_test.shape
    new_temp = [0]
    new_temp.extend(new_shape[1:])
    temp_data = np.ones(new_temp)
    for i in range(tensor_comp):
        model = models[i]
        strain_pred = model.predict(ms_list)
        for j in range(num_ms):
            pass
            # draw_strains_compare(strain_list_test[i,j,center],strain_pred[j,center],label=labels[i])
            # draw_differences([strain_list_test[i,j,center]-strain_pred[j,center]],['FE - MKS for Microstructure %d' % (j+1)])
        temp_data = np.append(temp_data, np.expand_dims(strain_pred,axis=0), axis=0)
    output_data = np.concatenate((output_data, temp_data), axis=1)
            
    f_output = open("data_for_analysis.csv", "w")
    elements = output_data[0,0]
    elements = elements.size
    old_shape = output_data.shape
    output_data = np.reshape(output_data, (old_shape[0], old_shape[1], elements))
    ## loop over elements
    for e in range(elements):
        ## loop over tensor values
        for i in range(old_shape[0]):
            ## loop over microstructures (first half are actual values, second half are predicted by MKS)
            for j in range(old_shape[1]):
                f_output.write("%e," % output_data[i,j,e])
        f_output.write("\n")
    f_output.close()
