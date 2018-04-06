#NOV 2016 Stefano Bosisio
#Script to write distres file for submitting host-guest simulations
#Usage: ~/sire.app/bin/python writedistres.py TOP CRD restraint.dat 
#restraint.dat is a file which specifies which atoms should be constrained:
#In restraint.dat write: C13,GLU53:CB = 5.0,10.0,1.0 (example)
#where C13 is the guest atom, GLU53 is the residue name and CB the atoms of the host
#5.0 is Req, 10.0 the K value, 1.0 D
#

import os,sys, random
import math
from Sire.Tools.OpenMMMD import *
from Sire.Tools import Parameter, resolveParameters



def createSystem(molecules):
    #print("Applying flexibility and zmatrix templates...")
    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system


def create_dictionary(restr_input):
    restr_dict = {}
    for lines in restr_input:
        if lines.startswith("#") or lines=="\n":
            continue
        else:
            splitter= lines.split("=") # ['MOL1:C5,GLY95:C5 ', ' 5.0,10.0,2.0']
            keys = splitter[0].split(",")
            lig = keys[0]
            host = keys[1]
            req = float(splitter[1].split(",")[0])
            keq = float(splitter[1].split(",")[1])
            D   = float(splitter[1].split(",")[2])
            dict_keys = "%s,%s" %(lig,host)
            restr_dict[dict_keys] = [req,keq,D]

    return restr_dict

def lighost_index(restr_dict,system):
    #restr_dict = ["MOL1:C1,GLY59:5"=[X,X,X]]
    #the ligand is usually the first residues
    lig_atoms = []
    ligand = system.molecule(MolNum(1)).atoms()
    for pairs in restr_dict:
        lig=pairs.split(",")[0]
        for at in ligand:
            if at.name().value()== lig:
                lig_atoms.append(at.index().value())
            else:
                continue

    host_atoms = []
    host_coords = []
    residues = system.molecule(MolNum(2)).residues()
    for pairs in restr_dict:
        res = "%s" % pairs.split(",")[1].split(":")[0]
        atom = "%s" %pairs.split(",")[1].split(":")[1]
        for residue in residues:
            res_name = residue.name().value()
            res_numb = int(residue.number().value())
            tofind = "%s%d"%(res_name,res_numb)
            if (tofind) == (res):
                print("FOUND")
                print(res)
                for at in residue.atoms():
                    if at.name().value() == atom.strip():
                        print("FOUND ATOM")
                        print(atom)
                        host_atoms.append(at.number().value())
                        coords = at.property("coordinates")
                        host_coords.append(coords)

                    else:
                        continue
            else:
                continue

    counter = 0
    output=open("distres","w")
    output.write("distance restraints dictionary = {")
    for pairs in restr_dict:
        if counter==len(restr_dict)-1:
            strings = "(%s,%s):(%s,%s,%s)}"%(lig_atoms[counter],host_atoms[counter],restr_dict[pairs][0],\
                                restr_dict[pairs][1],restr_dict[pairs][2])
            output.write(strings)
            restr_dict[pairs].append(host_coords[counter])
            counter+=1
        else:
            strings = "(%s,%s):(%s,%s,%s), "%(lig_atoms[counter],host_atoms[counter],restr_dict[pairs][0],\
                                        restr_dict[pairs][1],restr_dict[pairs][2])
            output.write(strings)

            restr_dict[pairs].append(host_coords[counter])
            counter+=1

    print(restr_dict)



    return(restr_dict,lig_atoms,host_atoms)

def defineIntegrationDomain(restr_dict):
    r"""Definition of the integration domain starting from host coordiantes
    Parameters
    ----------
    restr_dict : dictionary
                 restr_dict[lig_idx,host_idx] = ([req,D,K],[avgx,avgy,avgz])
                 where avgx,avgy and avgz are the average coordinates
    Returns
    ----------
    space : 3D array
            space = [(-x,-y,-z)(x,y,z)]
            where -x,-y and -z are the min values of the space and x,y and z the
            max values of the space
    """
    counter = 0
    for pairs in restr_dict:
        coords = restr_dict[pairs][3:]

        if counter==0:
            max_x = coords[0][0]
            min_x = coords[0][0]
            max_y = coords[0][1]
            min_y = coords[0][1]
            max_z = coords[0][2]
            min_z = coords[0][2]
            counter+=1
        else:
            if coords[0][0]> max_x:
                max_x = coords[0][0]
            if coords[0][0]<min_x:
                min_x = coords[0][0]
            if coords[0][1]> max_y:
                max_y = coords[0][1]
            if coords[0][1]< min_y :
                min_y = coords[0][1]
            if coords[0][2]>max_z:
                max_z = coords[0][2]
            if coords[0][2]<min_z :
                min_z = coords[0][2]

    #print(max_x,max_y,max_z,min_x,min_y,min_z)
    max_x += 5
    max_y += 5
    max_z += 5
    min_x -= 5
    min_y -= 5
    min_z -= 5
    #Adding a buffer region
    space = [(min_x,min_y, min_z), (max_x,max_y,max_z)]
    return space




#MAIN SCRIPT#

#take as input topology, coords and restraint.dat
topfile = sys.argv[1]
crdfile = sys.argv[2]
input_file = open(sys.argv[3],"r").readlines()
evaluation_bool = sys.argv[4] # if you want an evaluation of the restraint: True
                              # else just write the distres file
#create a restraint dictionary from restriant.dat
restr_dict = create_dictionary(input_file)
print("Atoms to be restrained...")
print(restr_dict)#sanity check

#Create a Sire system
amber = Amber()
(molecules, space) = amber.readCrdTop(crdfile, topfile)
system = createSystem(molecules)

#Extract the host atoms 
print("Wrinting distres file...")
restr_dict_2,lig_list,host_list  = lighost_index(restr_dict,system)

if evaluation_bool == "False":
    print("distres file written. Enjoy your simulation :)")
    sys.exit(-1)
else:
    print("Evaluation of the restraint...")
    #create the space of interation
    print("Creation of the domain of integration...")
    space=defineIntegrationDomain(restr_dict)
    #CONSTANTS
    delta = 0.10
    delta_over_two = delta/2.0
    deltavol = delta*delta*delta
    kb = 0.001987
    T = 298
    kbT = kb*T
    beta = 1/kbT
    ROT = 8 * pi**2
    Ztrans = 0.0
    Uavg = 0
    #
    #Grid creation
    Nx = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    Ny = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    Nz = int ( round ( ( space[1][0] - space[0][0] ) / delta ) )
    print("Number of elements to be evaluated %d" %(Nx*Ny*Nz))
    print("Evaluation...")
    for i in range(0,Nx):
        for j in range(0,Ny):
            for k in range(0,Nz):
                x = space[0][0] + delta*i + delta_over_two
                y = space[0][1] + delta*j + delta_over_two
                z = space[0][2] + delta*k + delta_over_two

                counter= 0
                for pairs in restr_dict:
                    #restr_dict[pairs]=[[req,K,D],[coords]]
                    x_dict = float(restr_dict[pairs][3][0])
                    y_dict = float(restr_dict[pairs][3][1])
                    z_dict = float(restr_dict[pairs][3][2])

                    distance = ((x-x_dict)**2 + (y-y_dict)**2 + (z-z_dict)**2)

                    upper_bound = (restr_dict[pairs][0]+ restr_dict[pairs][2])**2
                    intmd_bound = (restr_dict[pairs][0])**2
                    lower_bound = (restr_dict[pairs][0]- restr_dict[pairs][2])**2

                    if distance <= upper_bound and distance >= intmd_bound:
                        if counter == len(restr_dict)-1:
                            U = 0.0
                        else:
                            counter+=1

                    elif distance <= intmd_bound and distance >= lower_bound:
                        if counter== len(restr_dict)-1:
                            U = 0.0
                        else:
                            counter+=1
                    else:

                        dist = (math.sqrt(x**2 + y**2 + z**2))
                        K = (restr_dict[pairs][1])
                        D = (restr_dict[pairs][2])
                        U =(K*(dist-D)**2)


                        break

                Boltz = math.exp(-beta*U)*deltavol
                Ztrans += (Boltz)
                Uavg += U*Boltz*ROT
    #Calculation of Ztot, Uavg, S, Frestraint:
    Ztot = Ztrans*ROT
    Uavg /= (Ztot)

    Zideal = 1661.*ROT
    Delta_F = -kbT*math.log(Ztot/Zideal)
    minTDelta_S = -T*(kb*math.log(Ztot/Zideal)+Uavg/T)


    print ("Ztrans  = %8.5f Angstrom^3" % Ztrans)
    print ("Free energy Cost of removing the restraint = %8.5f kcal/mol" % -Delta_F)
    print("%8.5 kcal/mol is the value to be subtracted to the final free energy estimation" % Delta_F)

