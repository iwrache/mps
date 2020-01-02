import numpy as np
from numpy import linalg
Length=1

#左规范
def LeftCanonical(l,str,Dimension,flag,Matrix):
    ListA={}#用来临时存放A矩阵

    while flag!=-1:
        (Row,Column)=np.shape(Matrix)       
        U,S,V=linalg.svd(Matrix)
        SMatrix=np.zeros((Row,Column),dtype=np.float)
        #S在这里是一个向量，把它转化成矩阵形式
        S_Dimension=min(Row,Column)

        for i in range(min(Row,Column)):
            if S[i]==0:
                S_Dimension=i
                break
            SMatrix[i,i]=S[i]#将S转化为矩阵形式

        U=U[:,0:S_Dimension]
        SMatrix=SMatrix[0:S_Dimension,0:S_Dimension]
        V=V[0:S_Dimension,:]
        #如果S奇异值只有S_dimension个不为0，则截断U，S，V

        for j in range(Dimension):
            ListA[flag*Dimension+j]=U[j:j+Row-Dimension+1:Dimension,:]
               
        StimesV=np.mat(SMatrix)*np.mat(V)

        if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' ):
            (RowStimesV,ColumnStimesV)=np.shape(StimesV)
            Matrix=StimesV.reshape((Dimension*RowStimesV,-1))
            flag+=1

            if str=='mixed' and flag==l:
                return (ListA,SMatrix,V)
            #混合规范的情形

        else:
            print("extra number:%f"%StimesV)
            global LeftExtraNumber
            LeftExtraNumber=StimesV
            return ListA


#右规范
def RightCanonical(l,str,Dimension,flag,Matrix):
    ListB={}#临时存储B矩阵
    Length=flag
    while flag!=0:
        (Row,Column)=np.shape(Matrix)   
        U,S,V=linalg.svd(Matrix)
        SMatrix=np.zeros((Row,Column),dtype=np.float)
        S_Dimension=min(Row,Column)

        for i in range(min(Row,Column)):
            if S[i]==0:
                S_Dimension=i
                break
            SMatrix[i,i]=S[i]#将S转化为矩阵形式                             

        U=U[:,0:S_Dimension]
        SMatrix=SMatrix[0:S_Dimension,0:S_Dimension]
        V=V[0:S_Dimension,:]
        #如果S奇异值只有S_dimension个不为0，则截断U，S，V
        
        for j in range(Dimension):       
            ListB[(Length-flag)*Dimension+j]=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]
            #输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
            #与A矩阵输出稍有不同，A的跳着输出，B是连着输出

        UtimesS=np.matrix(U)*np.mat(SMatrix)

        if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
            (RowStimesV,ColumnStimesV)=np.shape(UtimesS)
            Matrix=UtimesS.reshape((-1,Dimension*ColumnStimesV))
            flag-=1
            
        elif str=='mixed' and flag==l:
            return (ListB,UtimesS)
        #混合规范的情形

        else:
            print("extra number:%d" %UtimesS)
            global RightExtraNumber
            RightExtraNumber=UtimesS
            return ListB
        

def MixedCanonical(l,Length,Dimension,Matrix):
    ListA=ListB={}#用来临时存放A和B
    (ListA,SMatrix,V)=LeftCanonical(l,'mixed',Dimension,0,Matrix)
    
    if l+2<=Length:
        NewMatrix=V.reshape((-1,Dimension))
        (ListB,UtimesS)=RightCanonical(l+2,'mixed',Dimension,Length,NewMatrix)
        
        (Row,Column)=np.shape(UtimesS)
        SiteLPlusOne=UtimesS.reshape((-1,Column*Dimension))
        (NewRow,NewColumn)=np.shape(SiteLPlusOne)
    else:
        SiteLPlusOne=V
        (NewRow,NewColumn)=np.shape(SiteLPlusOne)        
        #MixedCanonical没有多出来的数，因为它放到Smatrix里面了

    for j in range(Dimension):
        ListB[(Length-l-1)*Dimension+j]=SiteLPlusOne[:,j*int(NewColumn/Dimension):(j+1)*int(NewColumn/Dimension)]
        #输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
    return (ListA,SMatrix,ListB)


def OutputMPSMatrix(A,SMatrix,B,Dimension,Length,str,l):

    if str=='left':
        for i in range(Dimension*Length):
            print('第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(A[i])
        #print(A[2]*A[5]*A[8]*A[13]*A[19]*A[22]*LeftExtraNumber)
        #左数第n位第m个元素对应于A中的第(n-1)*Dimension+m-1个元素，如这里的对应于长度为6，维数为4的210132

    elif str=='right':
        for i in range(Dimension*Length):
            print('第%d位的第%d个矩阵：'%(Length-int(i/Dimension),i+1-Dimension*int(i/Dimension)))
            print(B[i])
        #print(B[22]*B[17]*B[12]*B[9]*B[7]*B[2]*RightExtraNumber)
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
        #将左规范和右规范结合起来即可，如这里对应的数是长度为6，维数为4,site为2的210132          


def Overlap(A,SMatrix,B,A_1,Dimension,Length,str,l):
    if str=='left':#之前是左规范的结果
        for i in range(Length):
            if i==0:
                BeforeResult=1#最里面一项的result是等于1的
            else:
                BeforeResult=result
            result=0#每次算完里面的都更新一下
            for j in range(Dimension):      
                result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*A[i*Dimension+j]#因为这里是实数，所以直接转置就行

    elif str=='right':
        for i in range(Length):
            if i==0:
                BeforeResult=1#最里面一项的result是等于1的
            else:
                BeforeResult=result
            result=0
            for j in range(Dimension):          
                result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*B[(Length-i-1)*Dimension+j]

    elif str=='mixed':
        for i in range(Length):
            if i==0:
                BeforeResult=1#最里面一项的result是等于1的
            else:
                BeforeResult=result
            result=0
            for j in range(Dimension):
                #因为mixed有个额外的S矩阵，所以多的那个extranumber就放到S矩阵里面了，所以相当于左规范和右规范都归一化了，所以向量和其本身内积是1，但mixed没有归一化，所以会多出来一个extra number
                if i<l:      
                    result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*A[i*Dimension+j]#因为这里是实数，所以直接转置就行
                elif i==l:

                    result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*SMatrix*B[(Length-i-1)*Dimension+j]  
                else:

                    result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*B[(Length-i-1)*Dimension+j]          
    
    print('overlap result:')
    print(result)
    print ('人家把overlap的结果算出来了的喵~你不夸夸人家吗~（害羞的低下头）')


def Reduced_Density_Matrix(A,SMatrix,B,Reduced_l,Dimension,Length):#还没有验证，不知道这么写对不对
    if str=='left':#之前是左规范的结果,因为左规范了，所以可以直接用书上121式
        for i in range(Length,Reduced_l,-1):
            if i==Length:
                BeforeResult=1
            else:
                BeforeResult=result
            result=0#每次算完里面的都更新一下
            for j in range(Dimension):      
                result+=A[(i-1)*Dimension+j]*BeforeResult*np.transpose(A[(i-1)*Dimension+j])#因为这里是实数，所以直接转置就行

    elif str=='right':#右规范了，所以直接用书上126式
        for i in range(Reduced_l):
            if i==0:
                BeforeResult=1#最里面一项的result是等于1的
            else:
                BeforeResult=result
            result=0
            for j in range(Dimension):          
                result+=np.transpose(B[(Length-i-1)*Dimension+j])*BeforeResult*B[(Length-i-1)*Dimension+j]

    elif str=='mixed':
        #这里套用127和128式，而这两个式子是一样的，因为中间是个对角阵
            result=SMatrix*np.mat(SMatrix)
            
    
    print('Reduced_Density_Matrix result:')
    print(result)
    print('')#暂时想不出来什么傲娇语录，先空着

def Generation_of_a_left_canonical_MPS(M,Length,Dimension):#左数第n位第m个矩阵存放在M[(n-1)*Dimension+m-1]中
    Left_A={}#用来存放左规范的mps
    for i in range(Length):
        if i==0:
            StimesV=1

        for j in range(Dimension):
            if j==0:
                New_U=StimesV*M[i*Dimension]
                
            else:
                New_U=np.hstack((New_U,StimesV*M[i*Dimension+j]))#将所有的M拼到一起，形成新的U，再将其奇异值分解

               
        New_U=New_U.reshape((Dimension*len(New_U),-1))

        (Row,Column)=np.shape(New_U)#相当于是左规范的一个逆过程，所以是先用hstack再reshape而不是vstack  
        S_Dimension=min(Row,Column)
        U,S,V=linalg.svd(New_U)
        SMatrix=np.zeros((Row,Column),dtype=np.float)#这里设置成局部变量就好

        for j in range(min(Row,Column)):
            if S[j]==0:
                S_Dimension=j
                break
            SMatrix[j,j]=S[j]#将S转化为矩阵形式
                                
        U=U[:,0:S_Dimension]
        SMatrix=SMatrix[0:S_Dimension,0:S_Dimension]
        V=V[0:S_Dimension,:]
        #如果S奇异值只有S_dimension个不为0，则截断U，S，V  
        StimesV=np.mat(SMatrix)*np.mat(V)

        for j in range(Dimension):
            Left_A[i*Dimension+j]=U[j:j+Row-Dimension+1:Dimension,:]     

    return Left_A


def Generation_of_a_right_canonical_MPS(M,Length,Dimension):#左数第n位第m个矩阵存放在M[(n-1)*Dimension+m-1]中
    Left_B={}#用来存放右规范的mps
    for i in range(Length-1,-1,-1):#Length-1到0
        if i==Length-1:
            UtimesS=1

        for j in range(Dimension):
            if j==0:
                New_U=M[i*Dimension]*UtimesS
                
            else:
                New_U=np.hstack((New_U,M[i*Dimension+j]*UtimesS))#将所有的M拼到一起，形成新的U，再将其奇异值分解

        
        (Row,Column)=np.shape(New_U)#相当于是右规范的一个逆过程
        U,S,V=linalg.svd(New_U)
        SMatrix=np.zeros((Row,Column),dtype=np.float)
        S_Dimension=min(Row,Column)

        for j in range(min(Row,Column)):
            if S[j]==0:
                S_Dimension=j
                break
            SMatrix[j,j]=S[j]#将S转化为矩阵形式
                                
        U=U[:,0:S_Dimension]
        SMatrix=SMatrix[0:S_Dimension,0:S_Dimension]
        V=V[0:S_Dimension,:]
        #如果S奇异值只有S_dimension个不为0，则截断U，S，V  
        UtimesS=np.matrix(U)*np.mat(SMatrix)

        for j in range(Dimension):
            Left_B[i*Dimension+j]=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]

    return Left_B


def SVD_compressing(List,Length,Dimension,Truncated_Dimension):#List为要压缩的矩阵列表,Truncated_Dimension表示截断的维度上限，Dimension是原来的维度
    M=Generation_of_a_left_canonical_MPS(List,Length,Dimension)#将要压缩的矩阵左规范，放入M中
    Left_B={}#用来存放右规范的mps
    for i in range(Length-1,-1,-1):#Length-1到0
        if i==Length-1:
            UtimesS=1

        for j in range(Dimension):
            if j==0:
                New_U=M[i*Dimension]*UtimesS
                
            else:
                New_U=np.hstack((New_U,M[i*Dimension+j]*UtimesS))#将所有的M拼到一起，形成新的U，再将其奇异值分解

        
        (Row,Column)=np.shape(New_U)#相当于是右规范的一个逆过程
        U,S,V=linalg.svd(New_U)
        SMatrix=np.zeros((Row,Column),dtype=np.float)
        S_Dimension=min(Row,Column)

        for j in range(min(Row,Column)):
            if S[j]==0:
                S_Dimension=j
                break
            SMatrix[j,j]=S[j]#将S转化为矩阵形式
                                
        U=U[:,0:S_Dimension]
        SMatrix=SMatrix[0:S_Dimension,0:S_Dimension]
        V=V[0:S_Dimension,:]
        #如果S奇异值只有S_dimension个不为0，则截断U，S，V  
  
        New_Dimension=min(Column,Row,Truncated_Dimension)
        U=U[:,0:New_Dimension]
        SMatrix=SMatrix[0:New_Dimension,0:New_Dimension]
        V=V[0:New_Dimension,:]
        #对应三种情况，分别是Trun<min(Col,Row)和Trun>min(Col,Row)&Col>Row和Trun>min(Col,Row)&Col<Row

        UtimesS=np.matrix(U)*np.mat(SMatrix)

        for j in range(Dimension):
            Left_B[i*Dimension+j]=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]

    return Left_B    


def Compress_a_matrix_product_state_iteratively(List,Length,Dimension,Truncated_Dimension):#List为要压缩的矩阵#应该有bug，还没调
    M=SVD_compressing(List,Length,Dimension,Truncated_Dimension)#SVD压缩后为一个右规范的MPS
    IterationTimes=1#设置最大迭代次数
    site=1#表明优化到了第site个

    while IterationTimes<=10:
        #if : #这里是达到了迭代条件，即书上150那个式子，我还没弄懂那个式子怎么推出来的，先空着
            #break
        if site==1:
            L=1
        if site==L:
            R=1
        
        for i in range(site):
            if i==0:
                BeforeResult=1#最里面一项的result是等于1的
            else:
                BeforeResult=L
            L=0#每次算完里面的都更新一下
            for j in range(Dimension):      
                L+=np.transpose(M[i*Dimension+j])*BeforeResult*List[i*Dimension+j]#因为这里是实数，所以直接转置就行
        
        for i in range(Length,site,-1):
            if i==0:
                BeforeResult=1#最里面一项的result是等于1的
            else:
                BeforeResult=L
            L=0#每次算完里面的都更新一下
            for j in range(Dimension):      
                R+=np.transpose(M[i*Dimension+j])*BeforeResult*List[i*Dimension+j]#因为这里是实数，所以直接转置就行    

        for i in range(Dimension):
            M[site*Dimension+i]=L*M[site*Dimension+i]*R

        site+=1
        IterationTimes+=1

    return M        
        

#以下是主函数部分
while(Length!=0):
    Length=int(input("Please input Length :"))
    if Length<=0:
        print("game over")
        break
    Dimension=int(input("Please input Dimension:"))
    if Dimension<=0:
        print("game over")
        break

    A=B={}#A表示左规范时的A矩阵,B表示右规范时的B矩阵

    str=input("left canonical or right canonical or mixed:")
    if str=='left': 
        Matrix=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
        SMatrix=0
        k=0
        print("please input matrix element")
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                #Matrix[i,j]=input()
                Matrix[i,j]=k
                k+=1          

        flag=0#flag+1表示第几位
        l=0
        A=LeftCanonical(0,'left',Dimension,flag,Matrix)#0在这里没有用，为mixed设的
        OutputMPSMatrix(A,0,0,Dimension,Length,'left',l)#3个0对应的位置是针对混合规范设的，这里没有用

    elif str=='right':
        Matrix=np.zeros((Dimension**(Length-1),Dimension),dtype=np.float)
        SMatrix=0
        k=0
        print("please input matrix element")
        for i in range(Dimension**(Length-1)):
            for j in range(Dimension):
                #Matrix[i,j]=input()
                Matrix[i,j]=k
                k+=1

        flag=Length#flag表示第flag个site
        l=0
        B=RightCanonical(0,'right',Dimension,flag,Matrix)#0在这里没有用，为mixed设的
        OutputMPSMatrix(0,0,B,Dimension,Length,'right',l)#3个0对应的位置是针对混合规范设的，这里没有用


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
            flag=0#flag+1表示第几位
            A=LeftCanonical(0,'left',Dimension,flag,Matrix)
            OutputMPSMatrix(A,0,0,Dimension,Length,'left',l)

        else:
            (A,SMatrix,B)=MixedCanonical(l,Length,Dimension,Matrix)
            OutputMPSMatrix(A,SMatrix,B,Dimension,Length,'mixed',l)

                
    else:
        print('哈？本小姐让你输入这个了吗？难得本小姐亲自为你服务你居然还这样，真是个不识趣的男人。')
        break

    str1=input('欧尼酱~人家想让你执行overlap，好不好嘛~（当然可以啊，我可爱的小公主/不好，滚）：')
    if str1=='当然可以啊，我可爱的小公主':
        print('那你得输入另外一个态前面的系数哦~')
        Matrix_1=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)

        k=0
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                #Matrix[i,j]=input()
                Matrix_1[i,j]=k
                k+=1

        A_1={}
        A_1=LeftCanonical(0,'left',Dimension,0,Matrix_1)
        Overlap(A,SMatrix,B,A_1,Dimension,Length,str,l)

    elif str1=='no':
        print('嘤嘤嘤，你又欺负人家（眼泪汪汪）')
    
    else:
        print('欧尼酱你个八嘎，居然输错了，人家明明都提醒过你了~(委屈巴巴)~人家不跟你玩了，哼~')
        break

    str2=input('Reduced density matrix(yes/no)：')
    if str2=='yes':
        Reduced_l=int(input('Please input reduced site l:'))
        if Reduced_l<=0 or Reduced_l>Length:
            print('error')
            break
        Reduced_Density_Matrix(A,SMatrix,B,Reduced_l,Dimension,Length)
    
    elif str2=='no':
        print('no')#此处傲娇语录先空着，回头一起补,以此类推
        #break

    else:
        print('哎呀，你个小笨蛋，都提示你怎么输入了你还能输错，再这样下去小心人家不陪你玩了哦。哼~（傲娇脸）')
        break

    str3=input('generation of a left canonical(yes/no):')
    if str3=='yes':
        Left_Canonical={}
        Left_Canonical=Generation_of_a_left_canonical_MPS(A,Length,Dimension)#拿已经左规范的A试一试
        for i in range(Dimension*Length):
            print('左规范后第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(Left_Canonical[i])
        #print(Left_Canonical[2]*Left_Canonical[5]*Left_Canonical[8]*Left_Canonical[13]*Left_Canonical[19]*Left_Canonical[22]*LeftExtraNumber)
        #左数第n位第m个元素对应于A中的第(n-1)*Dimension+m-1个元素，如这里的对应于长度为6，维数为4的210132

    elif str3=='no':
        print('no')

    else:
        break

    str4=input('generation of a right canonical(yes/no):')
    if str4=='yes':
        Right_Canonical={}
        Right_Canonical=Generation_of_a_right_canonical_MPS(A,Length,Dimension)
        for i in range(Dimension*Length):
            print('右规范后第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(Right_Canonical[i])
        #print(Right_Canonical[2]*Right_Canonical[5]*Right_Canonical[8]*Right_Canonical[13]*Right_Canonical[19]*Right_Canonical[22]*LeftExtraNumber)
        #左数第n位第m个元素对应于A中的第(n-1)*Dimension+m-1个元素，如这里的对应于长度为6，维数为4的210132
        #但为什么这里多一个负号？
    elif str4=='no':
        print('no')

    else:
        break

    str5=input('SVD compress(yes/no):')
    if str5=='yes':
        Truncated_Dimension=int(input('Please input upper bound:'))
        SVD_compressing_list=SVD_compressing(A,Length,Dimension,Truncated_Dimension)
        for i in range(Dimension*Length):
            print('SVD压缩后第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(SVD_compressing_list[i])
        #print(SVD_compressing_list[2]*SVD_compressing_list[5]*SVD_compressing_list[8]*SVD_compressing_list[13]*SVD_compressing_list[19]*SVD_compressing_list[22]*LeftExtraNumber)
        #左数第n位第m个元素对应于SVD_compressing_list中的第(n-1)*Dimension+m-1个元素，如这里的对应于长度为6，维数为4的210132

    elif str=='no':
        print('')

    else:
        break
