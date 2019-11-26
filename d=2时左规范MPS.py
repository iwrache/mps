import numpy as np
from numpy import linalg
Length=1

def svdfunction(flag,Matrix):
    (Row,Column)=np.shape(Matrix)
    
    U,S,V=linalg.svd(Matrix)
    SMatrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(min(Row,Column)):
            SMatrix[i,i]=S[i]
                            
    if Row>Column:
        U=U[:,0:Column]
        SMatrix=SMatrix[0:Column,:]
                        
    A0=U[0:int(Row/2),:]
    A1=U[int(Row/2):Row,:]
    print("A-%d-0"%(flag))
    print(A0)
    print("A-%d-1"%(flag))
    print(A1)          
    #print(np.mat(A0.T)*np.mat(A0)+np.mat(A1.T)*np.mat(A1))
  
    StimesV=np.mat(SMatrix)*np.mat(V)
    if min(Row,Column)!=1:
        MatrixBefore=StimesV[:,0:int(Column/2)]    
        MatrixAfter=StimesV[:,int(Column/2):Column]
        NewMatrix=np.vstack((MatrixBefore,MatrixAfter))
        #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
        flag+=1
        svdfunction(flag,NewMatrix)
    else:
        print("extra number:")
        print(StimesV)

        
        
while(Length!=0):
    Length=int(input("Please input Length :"))
    if Length<=0:
        print("game over")
        break  
    matrix=np.zeros((2,2**(Length-1)),dtype=np.float)
    k=1
    print("please input matrix element")
    for i in range(2):
       for j in range(2**(Length-1)):
           matrix[i,j]=float(k)
           k+=1
    flag=1#flag用于表明第几次执行svdfuncion
    svdfunction(flag,matrix)
