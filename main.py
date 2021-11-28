import torch
import numpy as np
import os
import torchani
import h5py
import random
import math
import sys
import parmed as pmd
from parmed.amber import Rst7

from functions import *

#----------------------------------------------------------------------------------------
#CARGO LISTA DE ESTRUCTURAS
#----------------------------------------------------------------------------------------
dicZ=[("H",1),("He",2),("Li",3),("Be",4),("B",5),("C",6),("N",7),("O",8),("F",9),("Ne",10)]
database = h5py.File('C:/Users/juans/OneDrive/Escritorio/Labo/ani1x-release.h5', 'r')
#database = h5py.File('/home/jota/BasesDeDatos/ani1x-release.h5', 'r')


os.remove("energias.txt")
file_object = open('energias.txt', 'a')
file_object.write('formula indice E_ref E_elec_cm5 E_elec_hir E_elec_mbis E_elec_mul_sol E_elec_mul_vac E_pol \n')
file_object.close()

lista = np.loadtxt('lista.txt', dtype=str)
lista = lista[254:]
nombre_archivos = []

for i in range(len(lista)): #iria len(lista)
    nombre_archivos.append(lista[i][0]+'_'+lista[i][1])

for numero, nombre_archivo in enumerate(nombre_archivos):
#----------------------------------------------------------------------------------------
#CARGO LOS .RST7 DE LAS ESTRUCTURAS SOLVATADAS
#----------------------------------------------------------------------------------------
   # r = parmed.load_file(nombre_archivo+'_sol_sp.rst7').coordinates[0]


#----------------------------------------------------------------------------------------
#EXTRAIGO LAS CARGAS DE LA BASE DE DATOS
#----------------------------------------------------------------------------------------
    species = database[lista[numero][0]+'/atomic_numbers']
    atoms_qm =  species.shape[0]
    formula = lista[numero][0]
    indice = int(lista[numero][1])

    q_cm5 = torch.Tensor(database[formula+'/wb97x_dz.cm5_charges'][indice])
    q_hir = torch.Tensor(database[formula+'/wb97x_dz.hirshfeld_charges'][indice])
    q_mbis = torch.Tensor(database[formula+'/wb97x_tz.mbis_charges'][indice])


#----------------------------------------------------------------------------------------
#EXTRAIGO LAS CARGAS MULLIKEN DE LIO
#----------------------------------------------------------------------------------------
    q_mul_sol = []
    inicio = 7 + atoms_qm
    fin = 8 + 2*atoms_qm
    with open(nombre_archivo+'_sol_MK', 'r') as f:
        for renglon, line in enumerate(f):
            if renglon > inicio and renglon < fin:
                q = line.split()[-1]
                q_mul_sol.append(float(q))
    q_mul_vac = []
    with open(nombre_archivo+'_vac_MK', 'r') as f:
        for renglon, line in enumerate(f):
            if renglon > inicio and renglon < fin:
                q = line.split()[-1]
                q_mul_vac.append(float(q))

#----------------------------------------------------------------------------------------
#ORDENO EL TENSOR DE COORDENADAS PARA QUE COINCIDA CON EL ORDEN DE LAS CARGAS
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
#ORDENO EL TENSOR DE CARGAS MULLIKEN
#----------------------------------------------------------------------------------------


    nombretop=nombre_archivo+"_vac_ord.prmtop" #revisar
    topfile = pmd.load_file(nombretop)
    atoms_qm=len(topfile.residues[0])
    listZ_top=[]
    listS_top=[]
    for i in range(0,atoms_qm):
        listZ_top.append(topfile.atoms[i].atomic_number)
        listS_top.append(dicZ[listZ_top[i]-1][0])


    listS_ord=sorted(listS_top)
    listZ_ord=[]
    for i in range(0,len(listS_ord)):
        for j in range(0,len(dicZ)):
            if (dicZ[j][0]==listS_ord[i]):
                listZ_ord.append(dicZ[j][1])

    a=''.join(listS_ord)

    q_mul_sol_ord = []
    q_mul_vac_ord = []

    rst7file=Rst7(nombre_archivo+'_sol_sp.rst7')

    coordenadas=rst7file.coordinates[:1][0]
    coordQM=coordenadas[0:atoms_qm]

    coordQM_ord=[]

    for i in range(0,len(listZ_ord)):
        j=0
        encontre=False
        while (j<len(listZ_top) and not encontre):
            if (listZ_top[j]==listZ_ord[i]):
                q_mul_sol_ord.append(q_mul_sol[j])
                q_mul_vac_ord.append(q_mul_vac[j])
                q_mul_sol = np.delete(q_mul_sol, j, axis=0)
                q_mul_vac = np.delete(q_mul_vac, j, axis=0)

                coordQM_ord.append(coordQM[j])
                coordQM = np.delete(coordQM, j,axis=0)

                listZ_top = np.delete(listZ_top, j,axis=0)
                encontre = True
            j=j+1

    # cargar en un tensor las coordenadas ordenadas
    coordQM_ord = torch.Tensor(np.array(coordQM_ord))

    coordQM_ord = torch.Tensor(np.array(coordQM_ord))
    coordMM = torch.Tensor(coordenadas[atoms_qm:])
    r = torch.cat((coordQM_ord, coordMM))


    q_mul_sol = torch.Tensor(q_mul_sol_ord)
    q_mul_vac = torch.Tensor(q_mul_vac_ord)


#----------------------------------------------------------------------------------------
#AGREGO LAS CARGAS DE LAS AGUAS
#----------------------------------------------------------------------------------------
# el número de aguas es un tercio de la diferencia entre la cantidad total de átomos y los del soluto
    wat = int((r.shape[0] -  atoms_qm)/3)
    q_cm5     = append_water_charges(q_cm5, wat)
    q_hir    = append_water_charges(q_hir, wat)
    q_mbis    = append_water_charges(q_mbis, wat)
    q_mul_vac = append_water_charges(q_mul_vac, wat)
    q_mul_sol = append_water_charges(q_mul_sol, wat)


#genero un array de dimension len(r) que tiene 0 o 1, segun el agua este o no en rango.
#para los atomosqm tiene un 0 pero eso no se usa nunca
    indice_a=get_indice_aguas_en_rango(r,atoms_qm,cutoff=15)

#----------------------------------------------------------------------------------------
#CALCULO LA ENERGIA QMMM DEL EMBEBIMIENTO
#----------------------------------------------------------------------------------------
    E_elec_cm5, E_elec_hir, E_elec_mbis, E_elec_mul_sol, E_elec_mul_vac, E_pol = compute_coulomb(r, q_cm5, q_hir, q_mbis, q_mul_sol, q_mul_vac, species, indice_a)

#----------------------------------------------------------------------------------------
#EXTRAIGO LA ENERGIA QMMM DE REFERENCIA
#----------------------------------------------------------------------------------------
    with open(nombre_archivo+'_sol_sp.out')as f:
        for line in f:
            if 'QM-MM elec.' in line:
                linea = line.split(" ")
                E_ref = float(linea[-1][:-2])*627.509391
#----------------------------------------------------------------------------------------
#GUARDO LAS ENERGIAS EN UN TXT
#----------------------------------------------------------------------------------------
    energias = str(E_ref)+' '+str(E_elec_cm5)+' '+str(E_elec_hir)+' '+str(E_elec_mbis)+' '+str(E_elec_mul_sol)+' '+str(E_elec_mul_vac)+' '+str(E_pol)
    file_object = open('energias.txt', 'a')
    file_object.write(formula+' '+str(indice)+' '+str(energias)+'\n')
    file_object.close()
