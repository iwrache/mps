
import numpy as np
from numpy import linalg
Row=1
while(Row!=0):
    Row=int(input("输入矩阵行数:"))
    if Row<=0:
        print("输入不合法，game over")
        break  
    Column=int(input("输入矩阵列数:"))
    if Column<=0:
        print("输入不合法，game over")
    print("输入矩阵元素：")
    Matrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(Row):
       for j in range(Column):
           Matrix[i,j]=float(input())
           
    (Row,Column)=np.shape(Matrix)
    U,S,V=linalg.svd(Matrix)
    SMatrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(min(Row,Column)):
            SMatrix[i,i]=S[i]
            
    if Row>Column:
        U=U[:,0:Column]
        SMatrix=SMatrix[0:Column,:]
    StimesV=np.mat(SMatrix)*np.mat(V)
    print("U")
    print(U)
    print("StimesV")
    print(StimesV)

'''
    NewMatrix=matrix.reshape((Row*2),-1)
    print("Matrix")
    print(matrix)
    print("NewMatrix")
    print(NewMatrix)
'''
'''
    U,S,V=linalg.svd(matrix)
    if Row>Column:
       U=U[:,0:Column]
    I=np.mat(U.T)*np.mat(U)
    print("I:")
    print(I)
'''
