import numpy as np
from numpy import linalg
Length=1

def svdfunction(Dimension,flag,Matrix):
    (Row,Column)=np.shape(Matrix)
    
    U,S,V=linalg.svd(Matrix)
    SMatrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(min(Row,Column)):
        SMatrix[i,i]=S[i]#将S转化为矩阵形式
                            
    if Row>Column:
        U=U[:,0:Column]
        SMatrix=SMatrix[0:Column,:]
        #如果原始矩阵的行大于列，只需取U的前Column列，S取前Column行即可
    for j in range(Dimension):
        print("第%d位的第%d个矩阵:"%(flag,j))
        A=U[j:j+Row-Dimension+1:Dimension,:]
        print(A)#输出A矩阵，也可以把它保存起来，这里没保存，直接输出了

        
    StimesV=np.mat(SMatrix)*np.mat(V)
    if min(Row,Column)!=1:
        (RowStimesV,ColumnStimesV)=np.shape(StimesV)
        NewMatrix=StimesV.reshape((Dimension*RowStimesV,-1))
        #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
        #调整一下输出格式就可以了，嗯，真香
        flag+=1
        svdfunction(Dimension,flag,NewMatrix)
    else:
        print("extra number:")
        print(StimesV)

        
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
    #k=0
    print("please input matrix element")
    for i in range(Dimension):
       for j in range(Dimension**(Length-1)):
           matrix[i,j]=input()
           #k+=1

    flag=1#flag用于表明第几次执行svdfuncion
    svdfunction(Dimension,flag,matrix)
#仅限L为偶数，L为奇数的时候好像有一点bug，还没调，调了再更新
