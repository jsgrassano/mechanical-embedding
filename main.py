#---------------------------------------------------------------------
import torch
import numpy as np
import h5py
import random
import math
import matplotlib.pyplot as plt
import parmed
from parmed.amber import Rst7

from functions import *

#---------------------------------------------------------------------
#CARGA LA BASE DE DATOS
database = h5py.File('C:/Users/juans/OneDrive/Escritorio/Labo/ani1x-release.h5', 'r')

# FORMULAS DE LOS COMPUESTOS QUE SE QUIEREN ANALIZAR
formulas = ["C2H5N1O2", "C3H7N1O2", "C5H11N1O2", "C6H13N1O2", "C6H13N1O2", "C9H11N1O2", "C9H11N1O3", "C11H12N2O2",
            "C3H7N1O3", "C4H9N1O3", "C6H9N3O2", "C6H14N2O2", "C6H14N4O2", "C4H8N2O3",  "C5H10N2O3",  "C5H9N1O2"]

#CANTIDAD DE ESTRUCTURAS A EXTRAER DE CADA FORMULA
N = 1000

#CANTIDAD DE AGUAS EN EL SISTEMA MM
wat = 4000

#---------------------------------------------------------------------


#---------------------------------------------------------------------
 
random.seed(40118)

#---------------------------------------------------------------------

energias = []
with open("energias.txt", "w") as energy_file:
    
    for formula in formulas:
        #EXTRAIGO N ESTRUCTURAS Y SUS 3 TIPOS DE CARGA
        index_list = random.sample(range(0,len(database[formula+'/coordinates'])),N)
        structures = extraer_N_estructuras(index_list, formula, database)
        charges_cm5, charges_hirshfeld, charges_mbis = extraer_N_cargas(index_list, formula, database)
        species = torch.Tensor(database[formula+'/atomic_numbers'])


    #AGREGO LAS CARGAS PARCIALES DE LAS AGUAS A CADA VECTOR
        charges_cm5_h20 = append_water_charges(charges_cm5, wat)
        charges_hirshfeld_h20 = append_water_charges(charges_hirshfeld, wat)
        charges_mbis_h20 = append_water_charges(charges_mbis, wat)
        
    #SOLVATO LAS ESTRUCTURAS
        solvated_structures = torch.empty(N, structures.shape[1] + 3*wat, 3)
        for i, structure in enumerate(structures):
            solvated_structures[i] = solvate(structure, wat)
            #GENERO EL ARCHIVO .rst7
            rst7 = Rst7(natom = len(solvated_structures[i]))
            rst7.coordinates = solvated_structures[i].numpy()
            fname=str(formula)+'_'+str(index_list[i])+'.rst7'
            rst7.write(fname)
    
    #DINÁMICA PARA OPTIMIZAR AGUAS
	

    #CALCULO LA ENERGÍA DE REFERENCIA
	
        E_ref


    #CALCULO LAS ENERGIAS DEL EMBEBIMIENTO MECÁNICO    
    
        for i in range(N):
            E_elec_cm5, E_pol_cm5 = compute_coulomb(solvated_structures[i], charges_cm5_h20[i], species)
            E_elec_hirshfeld, E_pol_hirshfeld = compute_coulomb(solvated_structures[i], charges_hirshfeld_h20[i], species)
            E_elec_mbis, E_pol_mbis = compute_coulomb(solvated_structures[i], charges_mbis_h20[i], species)
            print(formula, int(index_list[i]) , E_elec_cm5, E_pol_cm5,E_elec_hirshfeld, E_pol_hirshfeld,E_elec_mbis, E_pol_mbis, E_ref, file=energy_file)

energy_file.close()
