import numpy as np
from numpy import linalg
Length=1

#改成了While的形式，时间稍微少了一点，但数大了还是显示内存错误
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
        #这里应该是一个Dimension*Dimension的矩阵，应该右规范一下，还没有规范，先放着

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
    
    print('result:')
    print(result)

    print ('人家把overlap的结果算出来了的喵~你不夸夸人家吗~（害羞的低下头）')


    


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

            #print('都提示你怎么输入了还输入不对，怎么肥事啊小老弟。')

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
            #MixedWay=input('欧尼酱~您是想用递归算还是用while算的喵（>''<）~(While/递归):')

            (A,B)=MixedCanonical(l,Length,Dimension,Matrix)
            OutputMPSMatrix(Dimension,Length,'mixed',l)

                
    else:
        print('哈？本小姐让你输入这个了吗？难得本小姐亲自为你服务你居然还这样，真是个不识趣的男人。')
        break



    str1=input('欧尼酱~人家想让你执行overlap，好不好嘛~（当然可以啊，我可爱的小公主/不好，滚）：')
    if str1=='当然可以啊，我可爱的小公主':
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



    elif str1=='不好，滚':
        print('嘤嘤嘤，你又欺负人家（眼泪汪汪）')
    
    else:
        print('嘤嘤嘤，欧尼酱你个笨蛋，居然输错了，人家明明都提醒过你了~(委屈巴巴)~人家不跟你玩了，哼~')
        break



