import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from gmshProcMalha import *

# Parâmetros da geometria

Lx0 = -2
Lxf = 5
Ly0 = -0.5
Lyf = 0.5

dy= abs(Lyf - Ly0)

# Constantes da Função corrente

c1 = 0
c2 = 1

dc= abs(c2 - c1)

# Criação da malha

mshName= 'mshTri4.msh'

X = procMsh(mshName)[0]
Y = procMsh(mshName)[1]

npoints= len(X)

triang = mtri.Triangulation(X,Y)
IEN = procMsh(mshName)[2]
IENbound = procMsh(mshName)[3]
IENboundNames = procMsh(mshName)[4]
ne = IEN.shape[0]

# Definição das condições de contorno

Psicc = np.zeros( (npoints) , dtype='float')

cc1=[];cc2=[];cc3=[];cc4=[];cc5=[]
  
for elem in range(len(IENbound)):
  if IENboundNames[elem] == "inflow":
    cc1.append(IENbound[elem][0])
    cc1.append(IENbound[elem][1])
    Psicc[IENbound[elem][0]] = dc/dy*Y[IENbound[elem][0]] + c1 - dc/dy*Ly0
    Psicc[IENbound[elem][1]] = dc/dy*Y[IENbound[elem][1]] + c1 - dc/dy*Ly0
    
  if IENboundNames[elem] == "outflow":
    cc2.append(IENbound[elem][0])
    cc2.append(IENbound[elem][1])
    Psicc[IENbound[elem][0]] = dc/dy*Y[IENbound[elem][0]] + c1 - dc/dy*Ly0
    Psicc[IENbound[elem][1]] = dc/dy*Y[IENbound[elem][1]] + c1 - dc/dy*Ly0
    
  if IENboundNames[elem] == "airfWall":
    cc3.append(IENbound[elem][0])
    cc3.append(IENbound[elem][1])
    Psicc[IENbound[elem][0]] = (c1 + c2)/2
    Psicc[IENbound[elem][1]] = (c1 + c2)/2
    
  if IENboundNames[elem] == "bottomWall":
    cc4.append(IENbound[elem][0])
    cc4.append(IENbound[elem][1])
    Psicc[IENbound[elem][0]] = c1
    Psicc[IENbound[elem][1]] = c1
    
  if IENboundNames[elem] == "topWall":
    cc5.append(IENbound[elem][0])
    cc5.append(IENbound[elem][1])
    Psicc[IENbound[elem][0]] = c2
    Psicc[IENbound[elem][1]] = c2
  

cc = cc1+cc2+cc3+cc4+cc5

# Montagem das matrizes globais

K = np.zeros( (npoints,npoints),dtype='float' )
M = np.zeros( (npoints,npoints),dtype='float' )
Gx = np.zeros( (npoints,npoints),dtype='float' )
Gy = np.zeros( (npoints,npoints),dtype='float' )

# Loop de solução

for e in range(0,ne):
 v = IEN[e]

 # triangulo
 area = (1/2)* ( X[v[2]]*( Y[v[0]]-Y[v[1]]) \
               + X[v[0]]*( Y[v[1]]-Y[v[2]]) \
               + X[v[1]]*(-Y[v[0]]+Y[v[2]]) )

 bi = Y[v[1]]-Y[v[2]]
 bj = Y[v[2]]-Y[v[0]]
 bk = Y[v[0]]-Y[v[1]]

 ci = X[v[2]]-X[v[1]]
 cj = X[v[0]]-X[v[2]]
 ck = X[v[1]]-X[v[0]]

 kxele = (1.0/(4.0*area)) * np.array([ [bi*bi, bi*bj, bi*bk],
                                       [bj*bi, bj*bj, bj*bk],
                                       [bk*bi, bk*bj, bk*bk] ])
 kyele = (1.0/(4.0*area)) * np.array([ [ci*ci, ci*cj, ci*ck],
                                       [cj*ci, cj*cj, cj*ck],
                                       [ck*ci, ck*cj, ck*ck] ])
 kxyele = (1.0/(4.0*area)) * np.array([[bi*ci, bi*cj, bi*ck],
                                       [bj*ci, bj*cj, bj*ck],
                                       [bk*ci, bk*cj, bk*ck] ])
 
 gxele = (1.0/6.0) * np.array([ [bi, bj, bk],
                                [bi, bj, bk],
                                [bi, bj, bk] ])
 
 gyele = (1.0/6.0) * np.array([ [ci, cj, ck],
                                [ci, cj, ck],
                                [ci, cj, ck] ])

 kele = kxele + kyele

 mele = (area/12.0) * np.array([ [2.0, 1.0, 1.0],
                                  [1.0, 2.0, 1.0],
                                  [1.0, 1.0, 2.0] ])
 
 for ilocal in range(0,IEN.shape[1]):
  iglobal = IEN[e,ilocal]
  for jlocal in range(0,IEN.shape[1]):
   jglobal = IEN[e,jlocal]
   
   K[iglobal,jglobal] += kele[ilocal,jlocal]
   M[iglobal,jglobal] += mele[ilocal,jlocal]
   Gx[iglobal,jglobal] += gxele[ilocal,jlocal]
   Gy[iglobal,jglobal] += gyele[ilocal,jlocal]

# Imposição das condições de contorno

bPsi = np.zeros( (npoints) , dtype='float')

for i in cc:
 K[i,:] = 0 # zerando a linha i
 K[i,i] = 1.0 # 1.0 na diagonal
 bPsi[i]   = Psicc[i] # impondo T em b

# Solução da equação

Psi = np.linalg.solve(K,bPsi)

vx = np.linalg.solve(M,(Gy@Psi))
vy = np.linalg.solve(M,(-Gx@Psi))

v = np.sqrt(vx*vx + vy*vy)

# PLOTANDO OS RESULTADOS

#Malha
plt.triplot(X,Y,IEN,'ko-')
plt.plot(X[cc1],Y[cc1],'bo')
plt.plot(X[cc2],Y[cc2],'go')
plt.plot(X[cc3],Y[cc3],'yo')
plt.plot(X[cc4],Y[cc4],'ro')
plt.plot(X[cc5],Y[cc5],'mo')
plt.show()

#Resultados
fig1, (ax1,ax2,ax3,ax4) = plt.subplots(nrows = 4)

psiPlot = ax1.tricontour(triang, Psi , levels= 25 , cmap = 'viridis')
vxPlot = ax2.tricontourf(triang, vx , levels= 200 , cmap = 'viridis')
vyPlot = ax3.tricontourf(triang, vy , levels= 200 , cmap = 'viridis')
vPlot = ax4.tricontourf(triang, v , levels= 200 , cmap = 'viridis')

fig1.colorbar(psiPlot, ax=ax1)
fig1.colorbar(vxPlot, ax=ax2)
fig1.colorbar(vxPlot, ax=ax3)
fig1.colorbar(vPlot, ax=ax4)

for ax in (ax1,ax2,ax3,ax4):
  ax.plot(X[cc1],Y[cc1],'k-')
  ax.plot(X[cc2],Y[cc2],'k-')
  ax.plot(X[cc3],Y[cc3],'k-')
  ax.plot(X[cc4],Y[cc4],'k-')
  ax.plot(X[cc5],Y[cc5],'k-')

plt.show()