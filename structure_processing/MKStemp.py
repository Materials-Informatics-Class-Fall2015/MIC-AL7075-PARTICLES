import numpy as np
import os
import copy
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
import pymks.tools as pt
import ReadResponses as RR
import ReadMSFunction as RMS


def trainComponents(dir_train, load=1):
    ## train all components for the xx yy zz xy xz yz directions for this specific imposed strain
    strain_list_train = RR.readDirectory(dir_train)
    tensor_comp = strain_list_train.shape[0]
    
    n = strain_list_train.shape[3]
    center = (n - 1) / 2
    
    ms_delta = make_delta_microstructures(n_phases = 2, size = (n,n,n))
    # draw_microstructures(ms_delta[:,center])
    prim_basis = PrimitiveBasis(n_states = 2)
    
    # for strain in strain_list_train:
        # strain = strain[:,center,:,:]
        # pt.draw_strains(strain)
    
    models = []
    for i in range(tensor_comp):
        temp_train = strain_list_train[i]
        temp_train /= load
        model = MKSLocalizationModel(basis = prim_basis)
        model.fit(ms_delta,temp_train)
        coeff = model.coeff
        models.append(model)
        # draw_coeff(coeff[center])
    return models
    
def trainSingleLoadLevelModels(dir, start=0, load=.002):
    ## models returned are of the order macro load, tensor component, then the spatial model
    models = []
    for i in range(start, start+6):
        temp_dir = os.path.join(dir,str(i))
        temp_models = trainComponents(temp_dir, load=load)
        models.append(temp_models)
    ## reorder models now to match Voigt form (from 11, 12, 13, 22, 23, 33 to 11, 22, 33, 12, 13, 23)
    temp = models[1]
    models[1] = models[3]
    models[3] = temp
    temp = models[2]
    models[2] = models[5]
    temp2 = models[4]
    models[5] = temp2
    models[4] = temp
    return models
    
def predictArbitraryStrain(models, loads, ms_list):
    ## models ordered 11, 22, 33, 12, 13, 23 and normalized by the load, next level is tensor component
    ## loads are in same order
    
    n = ms_list.shape[1]
    center = (n - 1) / 2
    num_ms = ms_list.shape[0]
    tensor_comp = len(models[0])
    ## predicted values need to be in order of strain_component, MS, x, y, z
    results = np.zeros((tensor_comp,) + ms_list.shape)
    
    for i in range(tensor_comp):
        for j in range(len(loads)):
            model = models[j][i]
            temp_model = copy.deepcopy(model)
            temp_model.resize_coeff(ms_list[0].shape)
            strain_pred = temp_model.predict(ms_list)
            results[i,:] += strain_pred*loads[j]
            # print("Finished Predicting for large MS")
            # for j in range(num_ms):
                # temp = strain_pred[j,center]
                # temp = np.expand_dims(temp, axis=0)
                # print(temp.shape)
                # pt.draw_strains(temp)
    return results

if __name__ == "__main__":
    oldVersion = True


    if not oldVersion:
        import sys
        dir_train = sys.argv[-2]
        dir_test = sys.argv[-1]
        models = trainSingleLoadLevelModels(dir_train, load=-.002)
        ms_list = RMS.readDirectory(dir_test)
        loads = [-.00145,.004,-.00145,0,0,0]
        strain_pred = predictArbitraryStrain(models, loads, ms_list)
        strain_list_test = RR.readDirectory(dir_test)
        tensor_comp = strain_list_test.shape[0]
        num_ms = ms_list.shape[0]
        center = (ms_list.shape[1]-1) / 2
        labels = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
        for i in range(tensor_comp):
            for j in range(num_ms):
                draw_strains_compare(strain_list_test[i,j,center],strain_pred[i,j,center],label=labels[i])
        
    if oldVersion:
        import sys
        dir_test = os.getcwd()
        dir_predict = None
        dir_train = None
        if(os.path.isdir(sys.argv[-1])):
            dir_test = sys.argv[-1]
        dir_train = os.getcwd()
        if(os.path.isdir(sys.argv[-2])):
            dir_train = sys.argv[-2]
        if(len(sys.argv) > 2 and os.path.isdir(sys.argv[-3])):
            dir_predict = sys.argv[-1]
            dir_test = sys.argv[-2]
            dir_train = sys.argv[-3]
            
        labels = ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
        
        strain_list_train = RR.readDirectory(dir_train)
        tensor_comp = strain_list_train.shape[0]
        
        n = strain_list_train.shape[3]
        center = (n - 1) / 2
        
        ms_delta = make_delta_microstructures(n_phases = 2, size = (n,n,n))
        # draw_microstructures(ms_delta[:,center])
        prim_basis = PrimitiveBasis(n_states = 2)
        
        # for strain in strain_list_train:
            # strain = strain[:,center,:,:]
            # pt.draw_strains(strain)
        
        models = []
        for i in range(tensor_comp):
            temp_train = strain_list_train[i]
            model = MKSLocalizationModel(basis = prim_basis)
            model.fit(ms_delta,temp_train)
            coeff = model.coeff
            models.append(model)
            # draw_coeff(coeff[center])
            
        if(dir_predict):
            ms_list = RMS.readDirectory(dir_predict)
            print("Read in all large MS")
            n = ms_list.shape[1]
            center = (n - 1) / 2
            num_ms = ms_list.shape[0]
            for i in range(tensor_comp):
                model = models[i]
                temp_model = copy.deepcopy(model)
                temp_model.resize_coeff(ms_list[0].shape)
                strain_pred = temp_model.predict(ms_list)
                print("Finished Predicting for large MS")
                for j in range(num_ms):
                    temp = strain_pred[j,center]
                    temp = np.expand_dims(temp, axis=0)
                    print(temp.shape)
                    pt.draw_strains(temp)
        
        
        strain_list_test = RR.readDirectory(dir_test)
        print(strain_list_test.shape)
        num_ms = strain_list_test.shape[1]
        ms_list = RMS.readDirectory(dir_test)
        n = ms_list.shape[1]
        center = (n - 1) / 2
        output_data = strain_list_test
        new_shape = strain_list_test.shape
        new_temp = [0]
        new_temp.extend(new_shape[1:])
        temp_data = np.ones(new_temp)
        for i in range(tensor_comp):
            model = models[i]
            temp_model = copy.deepcopy(model)
            temp_model.resize_coeff(ms_list[0].shape)
            strain_pred = temp_model.predict(ms_list)
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
