import numpy as np
from numpy import linalg
Length=1

def LeftCanonical(C,l,str,Dimension,flag,Matrix):
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
        #print("第%d位的第%d个矩阵:"%(flag,j))
        A=U[j:j+Row-Dimension+1:Dimension,:]
        #print(A)#输出A矩阵，也可以把它保存起来，这里没保存，直接输出了
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

    #(a,b)=np.shape(A)
    #print("第%d位的矩阵有%d行%d列:"%(flag,a,b))
     
    StimesV=np.mat(SMatrix)*np.mat(V)
    if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
        (RowStimesV,ColumnStimesV)=np.shape(StimesV)
        NewMatrix=StimesV.reshape((Dimension*RowStimesV,-1))
        #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
        #调整一下输出格式就可以了，嗯，真香
        flag+=1
        return LeftCanonical(C,l,str,Dimension,flag,NewMatrix)
        #用递归的方法空间复杂度太高，一般用递归能实现的用for循环和栈也能实现，先不管这些，等我先把算法全部写完再慢慢优化（或者换个高级的电脑就不优化了）
        #用递归的话把局部变量设成全局变量会好一点

    elif str=='mixed' and flag==l:
        #print("S[%d]:" %l)
        #print(SMatrix)
        if Row<Column:
            V=V[0:Row,:]
            SMatrix=SMatrix[:,0:Row]
        return (np.mat(C)*np.mat(SMatrix),V)
        

    else:
        print("extra number:")
        print(StimesV)
        print("C")
        print(C*StimesV)
       

#改成了While的形式，时间稍微少了一点，但数大了还是显示内存错误
def LeftCanonical_While(l,str,Dimension,flag,Matrix):
    global A#A表示左规范时的A矩阵
    A={}
    global SMatrix
    while flag!=-1:
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
            #print("第%d位的第%d个矩阵:"%(flag,j))
            A[flag*Dimension+j]=U[j:j+Row-Dimension+1:Dimension,:]

        
        StimesV=np.mat(SMatrix)*np.mat(V)
        if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' ):
            (RowStimesV,ColumnStimesV)=np.shape(StimesV)
            Matrix=StimesV.reshape((Dimension*RowStimesV,-1))
            #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
            #调整一下输出格式就可以了，嗯，真香
            flag+=1



            if str=='mixed' and flag==l:
                #print("S[%d]:" %l)
                #print(SMatrix)
                if Row<Column:
                    V=V[0:Row,:]
                    SMatrix=SMatrix[:,0:Row]

                return V

        else:
            print("extra number:")
            print(StimesV)

            break


#右规范的While形式
def RightCanonical_While(l,str,Dimension,flag,Matrix):
    global B
    B={}
    Length=flag
    while flag!=0:
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
            #print("第%d位的第%d个矩阵:"%(flag,j))
            B[(Length-flag)*Dimension+j]=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
            #与A矩阵输出稍有不同，A的跳着输出，B是连着输出

        UtimesS=np.matrix(U)*np.mat(SMatrix)

        if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
            (RowStimesV,ColumnStimesV)=np.shape(UtimesS)
            Matrix=UtimesS.reshape((-1,Dimension*ColumnStimesV))
            flag-=1
            
        elif str=='mixed' and flag==l:
            return UtimesS

        else:
            print("extra number:%d" %UtimesS)
            break
        
def RightCanonical(C,l,str,Dimension,flag,Matrix):
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
        #print("第%d位的第%d个矩阵:"%(flag,j))
        B=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
        #print(B)#输出B矩阵，也可以把它保存起来，这里没保存，直接输出了
        #与A矩阵输出稍有不同，A的跳着输出，B是连着输出
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
 
    #(a,b)=np.shape(B)
    #print("第%d位的矩阵有%d行%d列:"%(flag,a,b))

    UtimesS=np.matrix(U)*np.mat(SMatrix)

    if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
        (RowStimesV,ColumnStimesV)=np.shape(UtimesS)
        NewMatrix=UtimesS.reshape((-1,Dimension*ColumnStimesV))
        flag-=1
        return RightCanonical(C,l,str,Dimension,flag,NewMatrix)

    elif str=='mixed' and flag==l:
        #print("C(B)")
        #print(C)
        return (C,UtimesS)

    else:
        print("extra number:%d" %UtimesS)
        print("C")
        print(C*UtimesS)



def MixedCanonical(l,Length,Dimension,Matrix):
    (C_A,V)=LeftCanonical(1,l,'mixed',Dimension,1,Matrix)
    #(a,b)=np.shape(C_A)
    #print("C_A row:%d C_A column:%d"%(a,b))
    if l+2<=Length:
        NewMatrix=V.reshape((-1,Dimension))
        (C_B,UtimesS)=RightCanonical(1,l+2,'mixed',Dimension,Length,NewMatrix)
        
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
        #print("第%d位的第%d个矩阵:"%(l+1,j))
        B=SiteLPlusOne[:,j*int(NewColumn/Dimension):(j+1)*int(NewColumn/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
        if j==0:
            #(a,b)=np.shape(B)
            #print("B row:%d B column:%d"%(a,b))
            print('result')
            if l+2<=Length:
                print(np.mat(C_A)*B*np.mat(C_B))
            else:
                print(np.mat(C_A)*B)
        #print(B)#输出B矩阵，也可以把它保存起来，这里没保存，直接输出了
        #与A矩阵输出稍有不同，A的跳着输出，B是连着输出

def MixedCanonical_While(l,Length,Dimension,Matrix):
    V=LeftCanonical_While(l,'mixed',Dimension,0,Matrix)

    if l+2<=Length:
        NewMatrix=V.reshape((-1,Dimension))
        UtimesS=RightCanonical_While(l+2,'mixed',Dimension,Length,NewMatrix)
        

        (Row,Column)=np.shape(UtimesS)
        SiteLPlusOne=UtimesS.reshape((-1,Column*Dimension))
        (NewRow,NewColumn)=np.shape(SiteLPlusOne)
    else:
        SiteLPlusOne=V
        (NewRow,NewColumn)=np.shape(SiteLPlusOne)        
        #这里应该是一个Dimension*Dimension的矩阵，应该右规范一下，还没有规范，先放着

    for j in range(Dimension):
        B[(Length-l-1)*Dimension+j]=SiteLPlusOne[:,j*int(NewColumn/Dimension):(j+1)*int(NewColumn/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
    




def OutputMPSMatrix(Dimension,Length,str):
    if str=='left':
        for i in range(Dimension*Length):
            print('第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(A[i])
        #print(A[2]*A[5]*A[8]*A[13]*A[19]*A[22])
        #左数第n位第m个元素对应于A中的第(n-1)*Dimension+m-1个元素，如这里的对应于长度为6，维数为4的210132
    elif str=='right':
        for i in range(Dimension*Length):
            print('第%d位的第%d个矩阵：'%(Length-int(i/Dimension),i+1-Dimension*int(i/Dimension)))
            print(B[i])
        #print(B[22]*B[17]*B[12]*B[9]*B[7]*B[2])
        #左数第n位第m个元素对应于A中的第(Length-n)*Dimension+m-1个元素，如这里的对应于长度为6，维数为4的210132
    elif str=='mixed':
        for i in range(l*Dimension):
            print('第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(A[i])

        print('中间的S矩阵为：')
        print(SMatrix)

        for i in range((Length-l)*Dimension):
            print('第%d位的第%d个矩阵：'%(Length-int(i/Dimension),i+1-Dimension*int(i/Dimension)))
            print(B[i])

        #print(A[2]*A[5]*SMatrix*B[12]*B[9]*B[7]*B[2])
        #将左规范和右规范结合起来即可，如这里对应的数是长度为6，维数为4的210132          



    


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
        k=0
        print("please input matrix element")
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                #Matrix[i,j]=input()
                Matrix[i,j]=k
                k+=1      
        
        LeftWays=input('你想用递归方法执行还是for循环执行？（递归/For）')
        if LeftWays=='递归':
            flag=1#flag用于表明第几次执行svdfuncion
            LeftCanonical(1,0,'left',Dimension,flag,Matrix)#第二个和第三个参数不重要，为mix设的
        elif LeftWays=='For':
            flag=0#flag+1表示第几位
            LeftCanonical_While(0,'left',Dimension,flag,Matrix)
            OutputMPSMatrix(Dimension,Length,'left')
        else :
            print('哎呀，你个小笨蛋，都提示你怎么输入了你还能输错，再这样下去小心人家不陪你玩了哦。哼~（傲娇脸）')

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
        RightWays=input('你想用递归方法执行还是for循环执行？（递归/While）')
        if RightWays=='递归':
            RightCanonical(1,0,'right',Dimension,flag,Matrix)#第二个和第三个参数不重要，为mix设的
        elif RightWays=='Whilw':
            RightCanonical_While(0,'right',Dimension,flag,Matrix)
            OutputMPSMatrix(Dimension,Length,'right')
        else:
            print('都提示你怎么输入了还输入不对，怎么肥事啊小老弟。')

    elif str=='mixed':
        Matrix=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
        k=0
        #print("please input matrix element")
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                #Matrix[i,j]=input()
                Matrix[i,j]=k
                k+=1      
        l=int(input("please input site l:"))
        if l<=0 or l>Length:
            print('game over')
        elif l==Length:
            print('FBI warning:This means executing left-canonical')
            flag=1#flag用于表明第几次执行svdfuncion
            LeftCanonical(1,0,'left',Dimension,flag,Matrix)#第二个和第三个参数不重要，为mix设的
        else:
            MixedWay=input('欧尼酱~您是想用递归算还是用while算的喵(While/递归):')
            if MixedWay=='递归':
                MixedCanonical(l,Length,Dimension,Matrix)
            elif MixedWay=='While':
                MixedCanonical_While(l,Length,Dimension,Matrix)
                OutputMPSMatrix(Dimension,Length,'mixed')
            else:
                print('嘤嘤嘤，欧尼酱你个笨蛋，居然输错了，人家明明都提醒过你了~')


    else:
        print("game over")