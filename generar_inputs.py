#---------------------------------------------------------------------
import h5py
import numpy as np
import subprocess
#---------------------------------------------------------------------

def ejecutar(comando):
    comando=comando.split()
    subprocess.call(comando)

#CARGA LA BASE DE DATOS
database = h5py.File('/home/santi/Descargas/ani1x-release.h5', 'r')
periodic_table = {'H':1,
                   'He':2,
                   'Li':3,
                   'Be':4,
                   'B':5,
                   'C':6,
                   'N':7,
                   'O':8,
                   'F':9,
                   'Ne':10,
                   'Na':11,
                   'Mg':12,
                   'Al':13,
                   'Si':14,
                   'P':15,
                   'S':16,
                   'Cl':17,
                   'Ar':18}

#with open("formulas.txt", "r") as f:
#    formulas = f.readlines()
path = '/home/santi/0-Pol/molecules/'
#GENERO LOS ARCHIVOS .XYZ DE TODAS LAS ESTRUCTURAS CON CADA FORMULA MOLECULAR
for formula in ['C3H7N1O3']: #in formulas  #ITERO SOBRE LAS FORMULAS MOLECULARES
    ejecutar('mkdir '+path+str(formula)) #GENERO UNA CARPETA PARA CADA FORMULA EN EL ARCHIVO formulas.txt  #GENERO UNA CARPETA PARA CADA FORMULA
    index_list = list(range(3)) # len(database[formula+'/coordinates']) #ITERO SOBRE TODAS LAS ESTRUCTURAS DE CADA FORMULA
    atoms =  database[formula+'/atomic_numbers'] #EXTRAIGO LOS NUMEROS ATÓMICOS
    for index in index_list: 
        structure = np.array(database[formula+'/coordinates'][index]) #GUARDO LA I-ESIMA ESTRUCTURA
        with open(path+formula+"/"+formula+"_"+str(index)+'.xyz', "w") as XYZ: #ABRO UN ARCHIVO VACIO
            XYZ.write(str(len(atoms))+'\n'+'\n')
            for i, atom in enumerate(atoms):
              atom_charac = list(periodic_table.keys())[list(periodic_table.values()).index(atom)] #EXTRAIGO EL SIMBOLO DEL ATOMO
              line = atom_charac+" "+str(float(structure[i][0]))+" "+str(float(structure[i][1]))+" "+str(float(structure[i][2]))+"\n"
              XYZ.write(line)
        XYZ.close()

        #generar un input de molcas
        with open(path+formula+"/"+formula+"_"+str(index)+'.input', "w") as MOLCAS_INPUT: #ABRO UN ARCHIVO VACIO
            MOLCAS_INPUT.write('&GATEWAY\n')
            MOLCAS_INPUT.write('coord='+path+formula+"/"+formula+"_"+str(index)+'.xyz'+'\n')
            MOLCAS_INPUT.write('basis=ANO-L\ngroup=NOSYM\n&SEWARD\n&SFC\nKSDFT=PBE0\nLOPROP')
        
        #correr molcas
        ejecutar('export MOLCAS=/home/santi/programas/molcas/build')
        ejecutar('export PATH=$PATH:$MOLCAS')
        ejecutar('pymolcas '+path+formula+"/"+formula+"_"+str(index)+'.input -f')
        
        #borrar input 
        ejecutar('rm '+path+formula+"/"+formula+"_"+str(index)+'.xyz')
        ejecutar('rm '+path+formula+"/"+formula+"_"+str(index)+'.input')
        


