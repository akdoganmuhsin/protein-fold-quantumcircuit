from operator import index
from tkinter import Label
from rdkit.Chem.rdmolfiles import MolFromPDBFile
import numpy as np
import math

num = 0
resDict = dict()
bounds = []

AAdict = dict()

def Extracter(file_name):
    matrix = []
    molecule = MolFromPDBFile(file_name, sanitize=False)
    conf = molecule.GetConformer()
    for i, atom in enumerate(molecule.GetAtoms()):
        positions = conf.GetAtomPosition(i)
        if not atom.GetMonomerInfo().GetResidueName() == "HOH":
            matrix.append([atom.GetMonomerInfo(), atom.GetSymbol(), positions.x, positions.y, positions.z])
    return matrix

def compute_distance(first, second):
    return math.sqrt((first[0]-second[0])**2+(first[1]-second[1])**2+(first[2]-second[2])**2)

def matrixToAminoAcid(matrix):
    result = []
    buffer = []
    num = 1
    for i in range(len(matrix)):
        row = matrix[i]
        if num == row[0].GetResidueNumber():
            buffer.append(row)
        else:
            try:
                result.append(atomToAminoAcid(buffer))
            except ZeroDivisionError:
                # print(buffer)
                pass
            buffer = []
            buffer.append(row)
            num = row[0].GetResidueNumber()
    return result

def atomToAminoAcid(buffer):
    x = 0
    y = 0
    z = 0
    for row in buffer:
        x += row[2]
        y += row[3]
        z += row[4]
    return [x/len(buffer), y/len(buffer), z/len(buffer), row[0]]

def splitMatrix(matrix):
    result = dict()
    chains = []
    for row in matrix:
        chain = row[0].GetChainId()
        if chain in chains:
            result[chain].append(row)
        else:
            chains.append(chain)
            result[chain] = []
            result[chain].append(row)
    return chains, result

def computeChain(matrix, control, theta):
    global num
    boundNum = 0
    aadict = dict()
    for i in range(len(matrix)):
        name = matrix[i][3].GetResidueName()
        if not name in AAdict.keys():
            AAdict[name] = 0
        AAdict[name] += 1
        if not name in aadict.keys():
            aadict[name] = 0
        aadict[name] += 1
        for j in range(i+theta, len(matrix)):
            distance = compute_distance(matrix[i], matrix[j])
            if distance < control:
                # print(matrix[i][3].GetResidueName(), matrix[j][3].GetResidueName())
                # print(distance)
                bound = "-".join(sorted([matrix[i][3].GetResidueName(), matrix[j][3].GetResidueName()]))
                if not bound in bounds:
                    bounds.append(bound)
                    resDict[bound] = 0
                resDict[bound] +=1
                num += 1
                boundNum += 1
    return aadict, boundNum

def getSequence(matrix):
    res = []
    for row in matrix:
        res.append(row[3].GetResidueName())
    return res

control = 5
theta = 50
matrix = Extracter("1jfl.pdb")
chains, matrixDict = splitMatrix(matrix)
for chain in chains:
    # print(getSequence(matrixToAminoAcid(matrixDict[chain])))
    computeChain(matrixToAminoAcid(matrixDict[chain]), control, theta)
    # print(num)

# print(resDict)
# print(bounds)
# print(AAdict)



# ACD=list(AAdict.keys())
# ACDArray = np.array(ACD)
# shape = ( 5, 4 )
# ACDArray.reshape( shape )

# print(ACDArray)



# ##Quantum>>

# from qiskit import QuantumCircuit, Aer, assemble
# qc = QuantumCircuit(5)
# svsim = Aer.get_backend('aer_simulator')
# qc.save_statevector()
# qobj = assemble(qc)
# final_state = svsim.run(qobj).result().get_statevector()
# # Print the statevector neatly:
# array_to_latex(final_state, prefix="\\text{Statevector = }")

