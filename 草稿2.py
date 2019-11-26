import numpy as np
from numpy import linalg
Length=1

def svdfunction(C,Dimension,flag,Matrix):
    (Row,Column)=np.shape(Matrix)
    
    U,S,V=linalg.svd(Matrix)
    SMatrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(min(Row,Column)):
            SMatrix[i,i]=S[i]
                            
    if Row>Column:
        U=U[:,0:Column]
        SMatrix=SMatrix[0:Column,:]
    for j in range(Dimension):
        #print("第%d位的第%d个矩阵:"%(flag,j))
        A=U[j:j+Row-Dimension+1:Dimension,:]
        #print(A)
      
        if flag==1 and j==2:
            C=C*A
  
        elif flag==2 and j==2:
            C=C*A
           
        elif flag==3 and j==1:
            C=C*A
        elif flag==4 and j==2:
            C=C*A
        elif flag==5 and j==0:
            C=C*A
        elif flag==6 and j==1:
            C=C*A
        elif flag==7 and j==1:
            C=C*A
        elif flag==8 and j==1:
            C=C*A
        elif flag==9 and j==2:
            C=C*A
        
    StimesV=np.mat(SMatrix)*np.mat(V)
    if min(Row,Column)!=1:
        (RowStimesV,ColumnStimesV)=np.shape(StimesV)
        NewMatrix=StimesV.reshape((Dimension*RowStimesV,-1))
        #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
        flag+=1
        svdfunction(C,Dimension,flag,NewMatrix)
    else:
        print("extra number:")
        print(StimesV)
        print("C")
        print(C*StimesV)
        
while(Length!=0):
    Length=int(input("Please input Length :"))
    if Length<=0:
        print("game over")
        break
    Dimension=int(input("Please input Demension:"))
    if Dimension<=0:
        print("game over")
        break
    matrix=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
    k=0
    #print("please input matrix element")
    for i in range(Dimension):
       for j in range(Dimension**(Length-1)):
           matrix[i,j]=k
           k+=1

    flag=1#flag用于表明第几次执行svdfuncion
    svdfunction(1,Dimension,flag,matrix)
