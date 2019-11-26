import numpy as np
from numpy import linalg
Length=1

def svdfunction(flag,Matrix):
    (RowMatrix,ColumnMatrix)=np.shape(Matrix)
    U,S,V=linalg.svd(Matrix)
    print("V1")
    print(np.mat(V.T)*np.mat(V))
    if(RowMatrix<ColumnMatrix):
        V=V[0:RowMatrix,:]
    
    print("V2")
    print(np.mat(V.T)*np.mat(V))
    NewRow=min(RowMatrix,ColumnMatrix)
    if NewRow!=1:
        SMatrix=np.zeros((NewRow,NewRow),dtype=np.float)
        for i in range(NewRow):
            SMatrix[i,i]=S[i]
        StimesV=np.mat(SMatrix)*np.mat(V)

        NewMatrix=StimesV.reshape((NewRow*2),-1)

        flag+=1
        svdfunction(flag,NewMatrix)

        
        
while(Length!=0):
    Length=int(input("Please input Length :"))
    if Length==0:
        print("game over")
        break  
    matrix=np.zeros((2,2**(Length-1)),dtype=np.float)
    print("please input matrix element")
    for i in range(2):
       for j in range(2**(Length-1)):
           matrix[i,j]=float(input())
    print(matrix)
    flag=1#flag用于表明第几次执行svdfuncion
    svdfunction(flag,matrix)
