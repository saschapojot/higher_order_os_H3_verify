import numpy as np
import pandas as pd
from datetime import datetime
import sys
from scipy.special import hermite
import copy
import pickle
from pathlib import Path
from scipy import sparse
# import json

#this script tests accuracy for S4c, H3

# python readCSV.py groupNum rowNum, then parse csv
if len(sys.argv)!=4:
    print("wrong number of arguments")


group=int(sys.argv[1])
rowNum=int(sys.argv[2])

testNum=int(sys.argv[3])
inPath="./combinedData/groupNew"+str(group)+"/row"+str(rowNum)+"/"
inParamFileName="./inParamsNew"+str(group)+".csv"
#read parameters from csv
dfstr=pd.read_csv(inParamFileName)
oneRow=dfstr.iloc[rowNum,:]

g0=float(oneRow.loc["g0"])
omegam=float(oneRow.loc["omegam"])
omegap=float(oneRow.loc["omegap"])
omegac=float(oneRow.loc["omegac"])
er=float(oneRow.loc["er"])#magnification
r=np.log(er)
thetaCoef=float(oneRow.loc["thetaCoef"])
theta=thetaCoef*np.pi
Deltam=omegam-omegap
e2r=er**2
lmd=(e2r-1/e2r)/(e2r+1/e2r)*Deltam



N2=80

L2=80
N1=2

L1=1/np.sqrt(2*omegac)*0.8

k2ValsAll=[]
for n2 in range(0,int(N2/2)):
    k2ValsAll.append(2*np.pi/(2*L2)*n2)
for n2 in range(int(N2/2),N2):
    k2ValsAll.append(2*np.pi/(2*L2)*(n2-N2))
k2ValsAll=np.array(k2ValsAll)

dx1=2*L1/N1
dx2=2*L2/N2
print("dx1="+str(dx1))
print("dx2="+str(dx2))

x1ValsAll=np.array([-L1+dx1*n1 for n1 in range(0,N1)])
x2ValsAll=np.array([-L2+dx2*n2 for n2 in range(0,N2)])
x1ValsAllSquared=x1ValsAll**2

dtEst = 1e-6
tFlushStart=0
tFlushStop=0.001
flushNum=4000

tTotPerFlush=tFlushStop-tFlushStart

stepsPerFlush=int(np.ceil(tTotPerFlush/dtEst))
dt=tTotPerFlush/stepsPerFlush

timeValsAll=[]
for fls in range(0,flushNum):
    startingInd = fls * stepsPerFlush
    for j in range(0,stepsPerFlush):
        timeValsAll.append(startingInd+j)

timeValsAll=np.array(timeValsAll)*dt

#to construct psi analytical
# den_vec=[]
matSpace=np.zeros((N1,N2),dtype=complex)
for n1 in range(0,N1):
    for n2 in range(0,N2):
        x1SqTmp = x1ValsAllSquared[n1]
        x2Tmp = x2ValsAll[n2]
        # den_vec.append(1/2*g0*np.sqrt(2/omegam)-g0*omegac*np.sqrt(2/omegam)*x1SqTmp)
        matSpace[n1,n2]=np.exp(1j*x2Tmp/(1/2*g0*np.sqrt(2/omegam)-g0*omegac*np.sqrt(2/omegam)*x1SqTmp))
# print((den_vec))
def psiAnalytical(t):
    psiTmp=matSpace*np.exp(-1j/omegap*(np.cos(omegap*t)-1))
    psiTmp /= np.linalg.norm(psiTmp, "fro")
    return psiTmp


outDir="./groupNew"+str(group)+"/row"+str(rowNum)+"/test"+str(testNum)+"S4c_H3Verify/"
Path(outDir).mkdir(parents=True, exist_ok=True)
#time-independent part of f1
f1_space_vec=[1/2*g0*np.sqrt(2/omegam)-g0*omegac*np.sqrt(2/omegam)*elem for elem in x1ValsAllSquared]
f1_space_vec=np.array(f1_space_vec)
evo_mat_space=np.array(np.outer(f1_space_vec,k2ValsAll),dtype=complex)#used to do evolution in k2 space

def evolution1Step(j,psi):
    """

    :param j: time step
    :param psi: wavefunction at the beginning of the time step j
    :return:
    """
    tj=timeValsAll[j]

    psi_tilde = np.fft.fft(psi, axis=1, norm="ortho")

    #evolution part1
    mat1=-1j*1/6*dt*evo_mat_space*np.sin(omegap*tj)
    mat1=np.exp(mat1)
    psi_tilde=mat1*psi_tilde

    #evolution part2
    mat2=-1j*2/3*dt*evo_mat_space*np.sin(omegap*(tj+1/2*dt))
    mat2=np.exp(mat2)
    psi_tilde=mat2*psi_tilde

    #evolution part3
    mat3=-1j*1/6*dt*evo_mat_space*np.sin(omegap*(tj+dt))
    mat3=np.exp(mat3)
    psi_tilde=mat3*psi_tilde


    psiNext= np.fft.ifft(psi_tilde, axis=1, norm="ortho")

    return psiNext


tEvoStart=datetime.now()

psiNumericalCurr=psiAnalytical(0)

# psiAnaCurr=psiAnalytical(0)

for fls in range(0,flushNum):
    tFlushStart = datetime.now()
    startingInd = fls * stepsPerFlush
    diffPerFlush=[]
    for st in range(0,stepsPerFlush):
        j = startingInd + st
        psiNumericalNext = evolution1Step(j, psiNumericalCurr)
        psiNumericalCurr = psiNumericalNext

        psiAnaCurr = psiAnalytical(timeValsAll[j] + dt)
        diffTmp = np.linalg.norm(psiNumericalCurr - psiAnaCurr, ord="fro")
        diffPerFlush.append(diffTmp)


    # outData = {"diff": diffPerFlush}
    outFile = outDir + "flush" + str(fls) + "diff.pkl"

    with open(outFile,"wb") as fptr:
        pickle.dump(diffPerFlush,fptr)

    tFlushEnd = datetime.now()
    print("flush "+str(fls)+" time: ",tFlushEnd-tFlushStart)


tEvoEnd=datetime.now()
print("evo time: ",tEvoEnd-tEvoStart)