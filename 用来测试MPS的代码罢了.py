import numpy as np
from numpy import linalg
Length=1

def LeftCanonical(C,Dimension,flag,Matrix):
    (Row,Column)=np.shape(Matrix)
    
    U,S,V=linalg.svd(Matrix)
    SMatrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(min(Row,Column)):
        SMatrix[i,i]=S[i]#将S转化为矩阵形式
                            
    if Row>Column:
        U=U[:,0:Column]
        SMatrix=SMatrix[0:Column,:]
        #如果原始矩阵的行大于列，只需取U的前Column列，S取前Column行即可
        #还没优化，S的奇异值为0的项不用算的，优化完了再更新
    for j in range(Dimension):
        print("第%d位的第%d个矩阵:"%(flag,j))
        A=U[j:j+Row-Dimension+1:Dimension,:]
        print(A)#输出A矩阵，也可以把它保存起来，这里没保存，直接输出了
        if flag==1 and j==1:
            C=C*A
        elif flag==2 and j==1:
            C=C*A
        elif flag==3 and j==1:
            C=C*A
        elif flag==4 and j==2:
            C=C*A
        elif flag==5 and j==2:
            C=C*A
        elif flag==6 and j==1:
            C=C*A
        elif flag==7 and j==1:
            C=C*A
        elif flag==8 and j==1:
            C=C*A
        elif flag==9 and j==2:
            C=C*A
        elif flag==10 and j==0:
            C=C*A
        
    StimesV=np.mat(SMatrix)*np.mat(V)
    if min(Row,Column)!=1:
        (RowStimesV,ColumnStimesV)=np.shape(StimesV)
        NewMatrix=StimesV.reshape((Dimension*RowStimesV,-1))
        #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
        #调整一下输出格式就可以了，嗯，真香
        flag+=1
        LeftCanonical(C,Dimension,flag,NewMatrix)
    else:
        print("extra number:")
        print(StimesV)
        print("C")
        print(C*StimesV)

def RightCanonical(C,Dimension,flag,Matrix):
    (Row,Column)=np.shape(Matrix)
    
    U,S,V=linalg.svd(Matrix)
    SMatrix=np.zeros((Row,Column),dtype=np.float)
    for i in range(min(Row,Column)):
        SMatrix[i,i]=S[i]#将S转化为矩阵形式
                            
    if Row<Column:
        V=V[0:Row,:]
        SMatrix=SMatrix[:,0:Row]
        #如果原始矩阵的行小于列，只需取V的前Row行，S取前Column列即可
        #还没优化，S的奇异值为0的项不用算的，优化完了再更新（S的奇异值是按大到小排列？）
    
    for j in range(Dimension):
        
        print("第%d位的第%d个矩阵:"%(flag,j))
        B=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
        print(B)#输出B矩阵，也可以把它保存起来，这里没保存，直接输出了
        if flag==1 and j==2:
            C=B*C         
        elif flag==2 and j==0:
            C=B*C          
        elif flag==3 and j==1:
            C=B*C         
        elif flag==4 and j==2:
            C=B*C           
        elif flag==5 and j==2:
            C=B*C 
        elif flag==6 and j==1:
            C=B*C
        elif flag==7 and j==1:
            C=B*C
        elif flag==8 and j==1:
            C=B*C
        elif flag==9 and j==2:
            C=B*C
        
    UtimesS=np.matrix(U)*np.mat(SMatrix)
    if min(Row,Column)!=1:
        (RowStimesV,ColumnStimesV)=np.shape(UtimesS)
        NewMatrix=UtimesS.reshape((-1,Dimension*ColumnStimesV))
        flag-=1
        RightCanonical(C,Dimension,flag,NewMatrix)
    else:
        print("extra number:%d" %UtimesS)
        print("C")
        print(C*UtimesS)


while(Length!=0):
    Length=int(input("Please input Length :"))
    if Length<=0:
        print("game over")
        break
    Dimension=int(input("Please input Dimension:"))
    if Dimension<=0:
        print("game over")
        break
    
    str=input("left canonical or right canonical:")
    if str=='left': 
        Matrix=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
        k=0
        print("please input matrix element")
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                #Matrix[i,j]=input()
                Matrix[i,j]=k
                k+=1
        
        flag=1#flag用于表明第几次执行svdfuncion
        LeftCanonical(1,Dimension,flag,Matrix)

    elif str=='right':
        Matrix=np.zeros((Dimension**(Length-1),Dimension),dtype=np.float)
        k=0
        print("please input matrix element")
        for i in range(Dimension**(Length-1)):
            for j in range(Dimension):
                #Matrix[i,j]=input()
                Matrix[i,j]=k
                k+=1
        flag=Length#flag表示第flag个site
        RightCanonical(1,Dimension,flag,Matrix)

    else:
        print("game over")
