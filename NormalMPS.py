import numpy as np
from numpy import linalg
Length=1


def LeftCanonical(l,str,Dimension,flag,Matrix):
    ListA={}#用来临时存放A矩阵
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
            ListA[flag*Dimension+j]=U[j:j+Row-Dimension+1:Dimension,:]
               
        StimesV=np.mat(SMatrix)*np.mat(V)

        if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' ):
            (RowStimesV,ColumnStimesV)=np.shape(StimesV)
            Matrix=StimesV.reshape((Dimension*RowStimesV,-1))
            #不能用reshape函数，因为reshape函数转化的格式和要求转化的格式不一样
            #调整一下输出格式就可以了，嗯，真香
            flag+=1

            if str=='mixed' and flag==l:
                if Row<Column:
                    V=V[0:Row,:]
                    SMatrix=SMatrix[:,0:Row]
                return (ListA,V)

        else:
            print("extra number:%f"%StimesV)#overlap时输入第二个矩阵也会输出这个，可以自己调一下
            global LeftExtraNumber
            LeftExtraNumber=StimesV
            return ListA


#右规范的While形式
def RightCanonical(l,str,Dimension,flag,Matrix):
    ListB={}#临时存储B矩阵
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
            ListB[(Length-flag)*Dimension+j]=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
            #与A矩阵输出稍有不同，A的跳着输出，B是连着输出

        UtimesS=np.matrix(U)*np.mat(SMatrix)

        if (min(Row,Column)!=1 and str!='mixed') or (str=='mixed' and flag!=l):
            (RowStimesV,ColumnStimesV)=np.shape(UtimesS)
            Matrix=UtimesS.reshape((-1,Dimension*ColumnStimesV))
            flag-=1
            
        elif str=='mixed' and flag==l:
            return (ListB,UtimesS)

        else:
            print("extra number:%d" %UtimesS)
            global RightExtraNumber
            RightExtraNumber=UtimesS
            return ListB
        

def MixedCanonical(l,Length,Dimension,Matrix):
    ListA=ListB={}#用来临时存放A和B
    (ListA,V)=LeftCanonical(l,'mixed',Dimension,0,Matrix)
    
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
        ListB[(Length-l-1)*Dimension+j]=SiteLPlusOne[:,j*int(NewColumn/Dimension):(j+1)*int(NewColumn/Dimension)]#输出的是j*int(Column/Dimension)到(j+1)*int(Column/Dimension)-1列
    return (ListA,ListB)


def OutputMPSMatrix(Dimension,Length,str,l):

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


def Overlap(Dimension,Length,str,l):
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

                    result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*PreviousSmatrix*B[(Length-i-1)*Dimension+j]  
                else:

                    result+=np.transpose(A_1[i*Dimension+j])*BeforeResult*B[(Length-i-1)*Dimension+j]          
    
    print('overlap result:')
    print(result)
    print ('人家把overlap的结果算出来了的喵~你不夸夸人家吗~（害羞的低下头）')


def Reduced_Density_Matrix(Reduced_l,Dimension,Length,Execute_overlap):
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
        if Execute_overlap==1:
            result=PreviousSmatrix*np.mat(PreviousSmatrix)#计算overlap的时候更新了一下Smatrix，所以这里用PreviousSmatrix  
        else:
            result=SMatrix*np.mat(SMatrix)#这个时候没有执行overlap，所以一开始输入的保存在SMatrix中               
    
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
        U,S,V=linalg.svd(New_U)
        SMatrix_Left=np.zeros((Row,Column),dtype=np.float)#这里设置成局部变量就好

        for j in range(min(Row,Column)):
            SMatrix_Left[j,j]=S[j]#将S转化为矩阵形式
                                
        if Row>Column:
            U=U[:,0:Column]
            SMatrix_Left=SMatrix_Left[0:Column,:]
            #如果原始矩阵的行大于列，只需取U的前Column列，S取前Column行即可
            #还没优化，S的奇异值为0的项不用算的，优化完了再更新    
        StimesV=np.mat(SMatrix_Left)*np.mat(V)

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
        SMatrix_Right=np.zeros((Row,Column),dtype=np.float)#这里设置成局部变量就好

        for j in range(min(Row,Column)):
            SMatrix_Right[j,j]=S[j]#将S转化为矩阵形式
                                
        if Row<Column:
            V=V[0:Row,:]
            SMatrix_Right=SMatrix_Right[:,0:Row]
            #如果原始矩阵的行小于列，只需取V的前Row行，S取前Column列即可
            #还没优化，S的奇异值为0的项不用算的，优化完了再更新（S的奇异值是按大到小排列？） 
        UtimesS=np.matrix(U)*np.mat(SMatrix_Right)

        for j in range(Dimension):
            Left_B[i*Dimension+j]=V[:,j*int(Column/Dimension):(j+1)*int(Column/Dimension)]

    return Left_B


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
    global A,B#A表示左规范时的A矩阵,B表示右规范时的B矩阵
    A=B={}

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

        flag=0#flag+1表示第几位
        l=0
        A=LeftCanonical(l,'left',Dimension,flag,Matrix)#l在这里没有用，为mixed设的
        OutputMPSMatrix(Dimension,Length,'left',l)
        #else :
            #print('哎呀，你个小笨蛋，都提示你怎么输入了你还能输错，再这样下去小心人家不陪你玩了哦。哼~（傲娇脸）')

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
        l=0
        B=RightCanonical(l,'right',Dimension,flag,Matrix)#l在这里没有用，为mixed设的
        OutputMPSMatrix(Dimension,Length,'right',l)


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
            OutputMPSMatrix(Dimension,Length,'left',l)

        else:
            (A,B)=MixedCanonical(l,Length,Dimension,Matrix)
            OutputMPSMatrix(Dimension,Length,'mixed',l)

                
    else:
        print('哈？本小姐让你输入这个了吗？难得本小姐亲自为你服务你居然还这样，真是个不识趣的男人。')
        break

    str1=input('欧尼酱~人家想让你执行overlap，好不好嘛~（当然可以啊，我可爱的小公主/不好，滚）：')
    if str1=='当然可以啊，我可爱的小公主':
        Execute_overlap=1#表明执行了overlap
        print('那你得输入另外一个矩阵的元素哦~')
        Matrix_1=np.zeros((Dimension,Dimension**(Length-1)),dtype=np.float)
        k=0
        for i in range(Dimension):
            for j in range(Dimension**(Length-1)):
                #Matrix[i,j]=input()
                Matrix_1[i,j]=k
                k+=1
        global A_1
        A_1={}
        if str=='mixed':
            global PreviousSmatrix
            PreviousSmatrix=SMatrix
        A_1=LeftCanonical(0,'left',Dimension,0,Matrix_1)
        Overlap(Dimension,Length,str,l)

    elif str1=='no':
        Execute_overlap=0#表明没有执行overlap
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
        Reduced_Density_Matrix(Reduced_l,Dimension,Length,Execute_overlap)
    
    elif str2=='no':
        print('no')#此处傲娇语录先空着，回头一起补,以此类推
        #break

    else:
        print('game over')
        break

    str3=input('generation of a left canonical(yes/no):')
    if str3=='yes':
        Left_Canonical={}
        Left_Canonical=Generation_of_a_left_canonical_MPS(A,Length,Dimension)#拿已经左规范的A试一试
        for i in range(Dimension*Length):
            print('左规范后第%d位的第%d个矩阵：'%(int(i/Dimension)+1,i+1-Dimension*int(i/Dimension)))
            print(Left_Canonical[i])
        print(Left_Canonical[2]*Left_Canonical[5]*Left_Canonical[8]*Left_Canonical[13]*Left_Canonical[19]*Left_Canonical[22]*LeftExtraNumber)
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
