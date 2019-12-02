import numpy as np
from numpy import linalg
Length=1

def LeftCanonical(l,str,Dimension,flag,Matrix):
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


    #(a,b)=np.shape(A)
    #print("第%d位的矩阵有%d行%d列:"%(flag,a,b))
     
    StimesV=np.mat(SMatrix)*np.mat(V)
    if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
        (RowStimesV,ColumnStimesV)=np.shape(StimesV)
        NewMatrix=StimesV.reshape((Dimension*RowStimesV,-1))
        #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
        #调整一下输出格式就可以了，嗯，真香
        flag+=1
        return LeftCanonical(l,str,Dimension,flag,NewMatrix)
        #用递归的方法空间复杂度太高，一般用递归能实现的用for循环和栈也能实现，先不管这些，等我先把算法全部写完再慢慢优化（或者换个高级的电脑就不优化了）
        #用递归的话把局部变量设成全局变量会好一点

    elif str=='mixed' and flag==l:
        print("S[%d]:" %l)
        print(SMatrix)
        if Row<Column:
            V=V[0:Row,:]
            SMatrix=SMatrix[:,0:Row]
        return V
        

    else:
        print("extra number:")
        print(StimesV)
        print("C")
        print(C*StimesV)
       

def RightCanonical(l,str,Dimension,flag,Matrix):
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
        #与A矩阵输出稍有不同，A的跳着输出，B是连着输出

 
    #(a,b)=np.shape(B)
    #print("第%d位的矩阵有%d行%d列:"%(flag,a,b))

    UtimesS=np.matrix(U)*np.mat(SMatrix)

    if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
        (RowStimesV,ColumnStimesV)=np.shape(UtimesS)
        NewMatrix=UtimesS.reshape((-1,Dimension*ColumnStimesV))
        flag-=1
        return RightCanonical(l,str,Dimension,flag,NewMatrix)

    elif str=='mixed' and flag==l:
        #print("C(B)")
        #print(C)
        return UtimesS

    else:
        print("extra number:%d" %UtimesS)
        print("C")
        print(C*UtimesS)

def MixedCanonical(l,Length,Dimension,Matrix):
    V=LeftCanonical(l,'mixed',Dimension,1,Matrix)
    #(a,b)=np.shape(C_A)
    #print("C_A row:%d C_A column:%d"%(a,b))
    if l+2<=Length:
        NewMatrix=V.reshape((-1,Dimension))
        UtimesS=RightCanonical(l+2,'mixed',Dimension,Length,NewMatrix)
        
        #(a,b)=np.shape(C_B)
        #print("C_B row:%d C_B column:%d"%(a,b))
        (Row,Column)=np.shape(UtimesS)
        SiteLPlusOne=UtimesS.reshape((-1,Column*Dimension))
        (NewRow,NewColumn)=np.shape(SiteLPlusOne)
    else:
        SiteLPlusOne=V
        (NewRow,NewColumn)=np.shape(SiteLPlusOne)        
        #这里应该是一个Dimension*Dimension的矩阵，应该右规范一下，还没有规范，先放着

    for j in range(Dimension):       
        print("第%d位的第%d个矩阵:"%(l+1,j))
        B=SiteLPlusOne[:,j*int(NewColumn/Dimension):(j+1)*int(NewColumn/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
        print(B)#输出B矩阵，也可以把它保存起来，这里没保存，直接输出了
        #与A矩阵输出稍有不同，A的跳着输出，B是连着输出
'''
    (a,b)=np.shape(B)
    print("第%d位的矩阵有%d行%d列:"%(l+1,a,b))
'''

while(Length!=0):
    Length=int(input("Please input Length :"))
    if Length<=0:
        print("game over")
        break
    Dimension=int(input("Please input Dimension:"))
    if Dimension<=0:
        print("game over")
        break
    
    str=input("left canonical or right canonical or mixed:")
    if str=='left': 
        Matrix=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
        #k=0
        print("please input matrix element")
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                Matrix[i,j]=input()
                #Matrix[i,j]=k
                #k+=1      
        flag=1#flag用于表明第几次执行svdfuncion
        LeftCanonical(0,'left',Dimension,flag,Matrix)#前两个参数不重要，为mix设的

    elif str=='right':
        Matrix=np.zeros((Dimension**(Length-1),Dimension),dtype=np.float)
        #k=0
        print("please input matrix element")
        for i in range(Dimension**(Length-1)):
            for j in range(Dimension):
                Matrix[i,j]=input()
                #Matrix[i,j]=k
                #k+=1
        flag=Length#flag表示第flag个site
        RightCanonical(0,'right',Dimension,flag,Matrix)#前两个参数不重要，为mix设的

    elif str=='mixed':
        Matrix=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
        #k=0
        #print("please input matrix element")
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                Matrix[i,j]=input()
                #Matrix[i,j]=k
                #k+=1      
        l=int(input("please input site l:"))
        if l<=0 or l>Length:
            print('game over')
        elif l==Length:
            print('FBI warning:This means executing left-canonical')
            flag=1#flag用于表明第几次执行svdfuncion
            LeftCanonical(0,'left',Dimension,flag,Matrix)
        else:
            MixedCanonical(l,Length,Dimension,Matrix)


    else:
        print("game over")