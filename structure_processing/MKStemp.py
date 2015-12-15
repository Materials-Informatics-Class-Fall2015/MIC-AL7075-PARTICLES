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
import cPickle


def writeOutputs(output_data):
    f_output = open("data_for_analysis.csv", "w")
    elements = output_data[0,0]
    elements = elements.size
    old_shape = output_data.shape
    output_data = np.reshape(output_data, (old_shape[0], old_shape[1], elements), order='F')
    ## loop over elements
    for e in range(elements):
        ## loop over tensor values
        for i in range(old_shape[0]):
            ## loop over microstructures (first half are actual values, second half are predicted by MKS)
            for j in range(old_shape[1]):
                f_output.write("%e," % output_data[i,j,e])
        f_output.write("\n")
    f_output.close()

def trainComponents(dir_train, load=1):
    ## train all components for the xx yy zz yz xz xy directions for this specific imposed strain
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
    models = [0]*6
    index = [0, 5, 4, 1, 3, 2]
    for i in range(start, start+6):
        temp_dir = os.path.join(dir,str(i))
        temp_models = trainComponents(temp_dir, load=load)
        # component = temp_models[index[i]]
        models[index[i]] = temp_models
    ## reorder models now to match Voigt form (from 11, 12, 13, 22, 23, 33 to 11, 22, 33, 23, 13, 12)
    return models
    
def predictArbitraryStrain(models, loads, ms_list):
    ## models ordered 11, 22, 33, 23, 13, 12 and normalized by the load, next level is tensor component
    ## loads are in same order
    
    n = ms_list.shape[1]
    center = (n - 1) / 2
    num_ms = ms_list.shape[0]
    tensor_comp = len(models)
    ## predicted values need to be in order of strain_component, MS, x, y, z
    results = np.zeros((tensor_comp,) + ms_list.shape)
    
    for i in range(tensor_comp):
        temp = models[i]
        if(type(temp)==type([])):
            for j in range(len(temp)):
                model = temp[j]
                temp_model = copy.deepcopy(model)
                temp_model.resize_coeff(ms_list[0].shape)
                strain_pred = temp_model.predict(ms_list)
                results[j,:] += strain_pred*loads[i]
        else:
            model = models[i]
            temp_model = copy.deepcopy(model)
            temp_model.resize_coeff(ms_list[0].shape)
            strain_pred = temp_model.predict(ms_list)
            results[i,:] = strain_pred*loads[i]
                # print("Finished Predicting for large MS")
                # for j in range(num_ms):
                    # temp = strain_pred[j,center]
                    # temp = np.expand_dims(temp, axis=0)
                    # print(temp.shape)
                    # pt.draw_strains(temp)
    return results

if __name__ == "__main__":

    oldVersion = False
    TrainPredict(oldVersion, sys.argv)
    
    
def TrainPredict(oldVersion, args, loads=[0,0,0,0,0,0.0035]):
    old_dir = os.getcwd()
    if not oldVersion:
        import sys
        dir_train = args[-2]
        dir_test = args[-1]
        models_f = "models.p"
        model_f = os.path.join(os.getcwd(), models_f)
        if(os.path.exists(model_f)):
            model_f = open(model_f, 'r')
            models = cPickle.load(model_f)
            model_f.close()
        else:
            models = trainSingleLoadLevelModels(dir_train, load=-.002)
            model_f = open(model_f, 'w')
            cPickle.dump(models, model_f)
            model_f.close()
        
        # loads = [0,0,0,.002,0,0]
        # loads = [-.00073*1.5,.002*1.5,-.00073*1.5,0,0,0]
        # loads = [-.000719,.00199,-.000718,4.45e-7,4.7e-8,-1.3e-7]
        # loads = [0,.00199,0,0,0,0]
        # loads = [.002*-.3595,.002,.002*-.3595,0,0,0]
        
        try:
            ms_list = RMS.readDirectory(dir_test)
            strain_list_test = RR.readDirectory(dir_test)
        except:
            strain_list_test = None
        if(strain_list_test is not None):
            strain_pred = predictArbitraryStrain(models, loads, ms_list)
            tensor_comp = strain_list_test.shape[0]
            num_ms = ms_list.shape[0]
            center = (ms_list.shape[1]-1) / 2
            labels = ['xx', 'yy', 'zz', 'yz', 'xz', 'xy']
            # for i in range(tensor_comp):
                # for j in range(num_ms):
                    # draw_strains_compare(strain_list_test[i,j,center],strain_pred[i,j,center],label=labels[i])
            output_data = np.concatenate((strain_list_test, strain_pred), axis=1)
            os.chdir(dir_test)
            #writeOutputs(output_data)
            os.chdir(old_dir)
            return output_data
            
        else:
            ms_list = RMS.readDirectory(dir_test)
            ms_list = ms_list[:,:100,:100,:100]
            strain_pred = predictArbitraryStrain(models, loads, ms_list)
            os.chdir(dir_test)
            #writeOutputs(strain_pred)
            os.chdir(old_dir)
            return strain_pred
            
    if oldVersion:
        import sys
        dir_test = os.getcwd()
        dir_predict = None
        dir_train = None
        if(os.path.isdir(args[-1])):
            dir_test = args[-1]
        dir_train = os.getcwd()
        if(os.path.isdir(args[-2])):
            dir_train = args[-2]
        if(len(args) > 2 and os.path.isdir(args[-3])):
            dir_predict = args[-1]
            dir_test = args[-2]
            dir_train = args[-3]
            
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
        os.chdir(dir_test)
        writeOutputs(output_data)   
        
