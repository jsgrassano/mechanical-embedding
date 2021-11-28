import numpy as np
import torch
import math
import random

periodic_table = {'H':1,
                   'HA':1,
                   'HB2':1,
                   'HB3':1,
                   'He':2,
                   'Li':3,
                   'Be':4,
                   'B':5,
                   'C':6,
                   'CA':6,
                   'CB':6,
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
                   'SG':16,
                   'Cl':17,
                   'Ar':18}

polarizabilities = {'H':4.50711,
                   'HA':4.50711,
                   'HB2':4.50711,
                   'HB3':4.50711,
                   'He':1.38375,
                   'Li':164.1125,
                   'Be':37.74,
                   'B':20.5,
                   'C':11.3,
                   'CA':11.3,
                   'CB':11.3,
                   'N':7.4,
                   'O':5.3,
                   'F':3.74,
                   'Ne':2.6611,
                   'Na':162.7,
                   'Mg':71.2,
                   'Al':57.8,
                   'Si':37.3,
                   'P':25,
                   'S':19.4,
                   'SG':19.4,
                   'Cl':14.6,
                   'Ar':11.083} #https://doi.org/10.1080/00268976.2018.1535143

def XYZtoTensor(XYZ):
    #abro el archivo
    with open(XYZ, 'r') as f:
        #cuento las lineas del archivo
        line_count = 0
        for i, line in enumerate(f):
            if i == 0:
            #me quedo con la cantidad de atomos y hago listas vacias de especies y coordenadas
                species = [0]*int(line)
                coordinates = [[0,0,0]]*int(line)
            line_count += 1
            if i in range(2, len(species)+2):
                species[i-2] = line.split()[0]
        f.close()

        frames = math.floor(line_count / (len(species) + 2))
        #genero una lista vacia con la cantidad de frames que hay
        tensor_coordinates = [0]*frames
        tensor_species = [0]*frames

        #transformo los atomos en sus indices de la tabla periodica
        for i, atom in enumerate(species):
            species[i] = periodic_table[atom]

        #lo convierto en un tensor de la dimension requerida
        species = torch.tensor(species).unsqueeze(0)

        #recorro las lineas de nuevo :s
    with open(XYZ, 'r') as f:
        for i, line in enumerate(f):
            for numero_de_frame, frame in enumerate(tensor):

                #genero una lista con los elementos de cada linea
                lista_linea = line.split()
                #guardo las coordenadas
                coordinates[i] = lista_linea[-3:]
                print(coordinates)
            frame = np.array(coordinates, dtype=float)

        #transformo las coordenadas a tensores de pytorch
        tensor = torch.tensor(tensor, requires_grad=True, dtype=torch.float32)

    return tensor_species, tensor_coordinates


def TensortoXYZ(species, tensor, nombre):
    XYZ = open (nombre,'w')
    num_atoms = len(species[0])
    frames = tensor.shape[0]
    symbols  = [] #vector de los simbolos de los atomos
    for i, atom in enumerate(species[0]):
        symbols.append(list(periodic_table.keys())[list(periodic_table.values()).index(atom)])
    for i in range(frames):
        XYZ.write(num_atoms+'\n'+'archivo generado con TensortoXYZ'+'\n')
        for j in range(num_atoms):
            line = symbols[j]+" "+str(float(tensor[i][j][0]))+" "+str(float(tensor[i][j][1]))+" "+str(float(tensor[i][j][2]))+"\n"
            XYZ.write(line)
    XYZ.close()
    return

def extraer_N_estructuras(lista_indices, formula, database):
    coordenadas = database[formula+'/coordinates']
    estructuras_seleccionadas = []
    for i in lista_indices:
        estructuras_seleccionadas.append(coordenadas[i])
    estructuras_sel = torch.tensor(estructuras_seleccionadas)
    return estructuras_sel

def extraer_N_cargas(indice, formula, database):

    cargas_cm5        = database[formula+'/wb97x_dz.cm5_charges'][indice]
    cargas_hirshfeld = database[formula+'/wb97x_dz.hirshfeld_charges'][indice]
    cargas_mbis       = database[formula+'/wb97x_tz.mbis_charges'][indice]

    return cargas_cm5, cargas_hirshfeld, cargas_mbis


def unit_vector(vector):
    """ Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """Finds angle between two vectors"""
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def x_rotation(vector,theta):
    """Rotates 3-D vector around x-axis"""
    R = np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0, np.sin(theta), np.cos(theta)]])
    return np.dot(R,vector)

def y_rotation(vector,theta):
    """Rotates 3-D vector around y-axis"""
    R = np.array([[np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta), 0, np.cos(theta)]])
    return np.dot(R,vector)

def z_rotation(vector,theta):
    """Rotates 3-D vector around z-axis"""
    R = np.array([[np.cos(theta), -np.sin(theta),0],[np.sin(theta), np.cos(theta),0],[0,0,1]])
    return np.dot(R,vector)

def makewaters(inflim,suplim,VdWdist):
    x = np.arange(inflim, suplim, VdWdist)
    waters=[]
    for i in range(len(x)):
        for j in range(len(x)):
            for k in range(len(x)):
                pos = (x[i],x[j],x[k])
                waters.append(pos)
    waters=np.array(waters)
    return waters

def mergezandcoordXYZ(symbols,coord,waters):
    print(len(symbols)+len(waters))
    print("Molecule solvated")
    for i in range(0,len(simbolos)):
        print(symbols[i],coord[i][0].item(),coord[i][1].item(),coord[i][2].item())
    for i in range(0,int(len(waters)/3)):
        print('O ',waters[3*i][0].item(),waters[3*i][1].item(),waters[3*i][2].item())
        print('H ',waters[(3*i)+1][0].item(),waters[(3*i)+1][1].item(),waters[(3*i)+1][2].item())
        print('H ',waters[(3*i)+2][0].item(),waters[(3*i)+2][1].item(),waters[(3*i)+2][2].item())

def deleteOverlaped(waters,coord,VdW):
    newwaters=[]
    for i in range(0,len(waters)):
        tooclose=False
        for j in range(0,len(coord)):
          x=coord[j][0]-waters[i][0]
          y=coord[j][1]-waters[i][1]
          z=coord[j][2]-waters[i][2]
          dist=math.sqrt(x**2+y**2+z**2)
          if (dist < VdW):
            tooclose=True
        if not tooclose:
            newwaters.append(waters[i])
    newwaters=np.array(newwaters)
    return newwaters

def keepOnlyNCloser(waters,coord,N):
    if N >= len(waters):
        return waters
    watersWindex=[]
    for i in range(0,len(waters)):
        x0=coord[0][0]-waters[0][0]
        y0=coord[0][1]-waters[0][1]
        z0=coord[0][2]-waters[0][2]
        distMIN=math.sqrt(x0**2+y0**2+z0**2)
        for j in range(0,len(coord)):
            x=coord[j][0]-waters[i][0]
            y=coord[j][1]-waters[i][1]
            z=coord[j][2]-waters[i][2]
            dist=math.sqrt(x**2+y**2+z**2)
            if dist<distMIN:
               distMIN=dist
        vec=[waters[i][0],waters[i][1],waters[i][2],distMIN]
        watersWindex.append(vec)
    watersWindex=np.array(watersWindex)
    watersWindex=watersWindex[watersWindex[:, 3].argsort()]
    newwaters=[]
    for i in range(0,N):
        vec=(watersWindex[i][0],watersWindex[i][1],watersWindex[i][2])
        newwaters.append(vec)
    newwaters=np.array(newwaters)
    return newwaters

def randnormalvec():
    vec=np.array([(0.5-random.random()),(0.5-random.random()),(0.5-random.random())])
    norm=math.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    return vec/norm

def addHydrogens(waters):
    theta_WAT=109
    newwaters=[]
    for i in range(0,len(waters)):
        r=0.94
        vecH1=r*randnormalvec()
        vecH2=r*randnormalvec()

        x1, y1, z1 = vecH1
        r1=math.sqrt(x1**2+y1**2+z1**2)
        phi1=math.atan(y1/x1)
        theta1=math.acos(z1/r1)

        x2, y2, z2 = vecH2
        r2=math.sqrt(x2**2+y2**2+z2**2)
        phi2=math.atan(y2/x2)
        theta2=math.acos(z2/r2)
        #print(np.degrees(theta2))
        # print("ANTES",np.degrees(angle_between(vecH1,vecH2)))
        if (x1>0):
            rotx=-phi1
            roty=-theta1
        else:
            rotx=-phi1+np.radians(180)
            roty=-theta1

        vecH1=z_rotation(vecH1,rotx)
        vecH2=z_rotation(vecH2,rotx)
        vecH1=y_rotation(vecH1,roty)
        vecH2=y_rotation(vecH2,roty)

        theta2=np.radians(theta_WAT)
        x2=r2*math.sin(theta2)*math.cos(phi2)
        y2=r2*math.sin(theta2)*math.sin(phi2)
        z2=r2*math.cos(theta2)
        vecH2 = x2, y2, z2

        vecH1=z_rotation(vecH1,-rotx)
        vecH2=z_rotation(vecH2,-rotx)
        vecH1=y_rotation(vecH1,-roty)
        vecH2=y_rotation(vecH2,-roty)

        vecH1=vecH1+waters[i]
        vecH2=vecH2+waters[i]

        # check=True
        # if (abs(np.degrees(angle_between(vecH1-waters[i],vecH2-waters[i]))-109)>0.01):
        #     check=check and False

        newwaters.append(waters[i])
        newwaters.append(vecH1)
        newwaters.append(vecH2)
    newwaters=np.array(newwaters)

    return newwaters

def solvate(structure, wat):
    waters = makewaters(-30,30,3.1)
    waters = deleteOverlaped(waters, structure, 4)
    waters = keepOnlyNCloser(waters, structure, wat)
    waters = addHydrogens(waters)
    waters = torch.Tensor(waters)
    solvated_structure = torch.stack((structure, waters))

    return solvated_structure

def append_water_charges(charges, wat):
    charges = charges.numpy()
    charges_h20 = np.array([-0.834,0.417,0.417]*wat)
    total_charges = np.concatenate((charges, charges_h20))
    total_charges= torch.Tensor(total_charges)
    return total_charges


def dist(x1,y1,z1,x2,y2,z2):
    import math
    dist_cuad = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    dist = math.sqrt(dist_cuad)
    return dist

def E_field2(r,q,nqm,ntot,i):
    E_field=torch.zeros(3)
    for j in range(nqm,ntot):
        rij = np.array([r[j][0].item() - r[i][0].item(), r[j][1].item() - r[i][1].item(), r[j][2].item() - r[i][2].item()])
        distij = dist(r[i][0].item(), r[i][1].item(),r[i][2].item(), r[j][0].item(), r[j][1].item(),r[j][2].item())
        rij = rij/distij
        for k in range(3):
            E_field[k] = E_field[k]+(q[j].item()/distij**2)*rij[k]
    E_field2=E_field[0]**2+E_field[1]**2+E_field[2]**2
    return E_field2

def estaEnRango(r,nqm,j,cutoff):
    esta=False
    i=0
    while (not esta and i<nqm):
        rij = np.array([r[j][0].item() - r[i][0].item(), r[j][1].item() - r[i][1].item(), r[j][2].item() - r[i][2].item()])
        distij = dist(r[i][0].item(), r[i][1].item(),r[i][2].item(), r[j][0].item(), r[j][1].item(),r[j][2].item())
        esta = esta or distij < cutoff
        i=i+1
    return esta

def get_indice_aguas_en_rango(r,atoms_qm,cutoff):
    lista=[]
    for i in range(0,atoms_qm):
        lista.append(0)
    for i in range(atoms_qm,len(r)):
        if (estaEnRango(r,atoms_qm,i,cutoff)):
            lista.append(1)
        else:
            lista.append(0)
    return lista


def compute_coulomb(r, q_cm5, q_hir, q_mbis, q_mul_sol, q_mul_vac, species, lista):
    nqm = len(species)
    ntot = len(q_cm5)
    E_elec_cm5 = 0.0
    E_elec_hir = 0.0
    E_elec_mbis = 0.0
    E_elec_mul_sol = 0.0
    E_elec_mul_vac = 0.0
    E_pol = 0.0
    for i in range(nqm):
        qi_cm5 = q_cm5[i].item()
        qi_hir = q_hir[i].item()
        qi_mbis = q_mbis[i].item()
        qi_mul_sol = q_mul_sol[i].item()
        qi_mul_vac = q_mul_vac[i].item()
        alpha = 0.14818471*polarizabilities[list(periodic_table.keys())[list(periodic_table.values()).index(species[i])]]
        E_field=torch.zeros(3)
        for j in range(nqm, ntot):
            qj  = q_cm5[j].item()
            rij = np.array([r[j][0].item() - r[i][0].item(), r[j][1].item() - r[i][1].item(), r[j][2].item() - r[i][2].item()])
            distij = dist(r[i][0].item(), r[i][1].item(),r[i][2].item(), r[j][0].item(), r[j][1].item(),r[j][2].item())
            rij = rij/distij
            #if estaEnRango(r,nqm,j,cutoff):
            if lista[j]==1:
                E_elec_cm5 = E_elec_cm5 + qi_cm5*qj/distij
                E_elec_hir = E_elec_hir + qi_hir*qj/distij
                E_elec_mbis = E_elec_mbis + qi_mbis*qj/distij
                E_elec_mul_sol = E_elec_mul_sol + qi_mul_sol*qj/distij
                E_elec_mul_vac = E_elec_mul_vac + qi_mul_vac*qj/distij
                for k in range(3):
                    E_field[k] = E_field[k]+(qj/distij**2)*rij[k]
        E_field2 = E_field[0]**2+E_field[1]**2+E_field[2]**2
        E_pol = E_pol + alpha*E_field2
    E_elec_cm5 = E_elec_cm5*0.529177249*627.509391
    E_elec_hir = E_elec_hir*0.529177249*627.509391
    E_elec_mbis = E_elec_mbis*0.529177249*627.509391
    E_elec_mul_sol = E_elec_mul_sol*0.529177249*627.509391
    E_elec_mul_vac = E_elec_mul_vac*0.529177249*627.509391
    E_pol = -0.5*E_pol.item()*0.529177249*627.509391

    return E_elec_cm5, E_elec_hir, E_elec_mbis, E_elec_mul_sol, E_elec_mul_vac, E_pol
