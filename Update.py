import numpy as np
from numpy import ones,einsum
from scipy.linalg import norm, svd
import math
import pandas as pd

class Update:
    def __init__(self,tensor,center,Number,DataClass,Dimension,bond_dimension,maxbond,cutoff = 0.01):

        self.tensor=tensor#数据点的mps
        self.center=center#中心点的mps
        self.Number=Number#一共有Number个数
        self.DataClass=DataClass#将数据分为DataClass个类别
        self.Dimension=Dimension#数据一共有多少维
        self.current_bond = 0#当前更新的bond编号
        self.cumulants = {}#存储
        self.NumberClass = {}#数据类别
        self.bond_dimension=[bond_dimension for i in range(self.DataClass)]#中心点之间bond维数
        self.cutoff=0.01#截断比例
        self.merged_matrix=[None for i in range(self.DataClass)]#合并起来的4阶张量
        self.maxbond=maxbond#最大的bond维数

    def right_cano(self):
        for i in range(self.DataClass):
            for j in range(self.Dimension-1,-1,-1):
                dl = self.center[i][j].shape[0]
                dr = self.center[i][j].shape[2]
                U,s,V=svd(self.center[i][j].reshape(dl,dr*2).T)
                myd = min(dl,dr*2)
                Q=U[:,:myd]
                #R=np.diag(s)@(V[:,:myd].T)
                self.center[i][j]=Q.T.reshape(-1,2,dr)
                #self.center[i][j-1]=np.einsum("abc,ci->abi",self.center[i][j-1],R.T)

        #self.judge_left_cano()
        #self.judge_right_cano()

    def judge_left_cano(self,Dir):
        if Dir =='r':
            for j in range(self.DataClass):
                a = np.ones((1, 1))
                for i in range(self.current_bond):
                    a = np.einsum('mi,ijk,njm->nk', a, self.center[j][i], self.center[j][i].T)
                if (a-np.identity(a.shape[0])>0.001).any():
                    print('left_a:', a)
        else:
            for j in range(self.DataClass):
                a = np.ones((1, 1))
                for i in range(self.current_bond+1):
                    a = np.einsum('mi,ijk,njm->nk', a, self.center[j][i], self.center[j][i].T)
                if (a-np.identity(a.shape[0])>0.001).any():
                    print('left_a:', a)

    def judge_right_cano(self,Dir):
        if Dir =='r':
            for j in range(self.DataClass):
                a = np.ones((1, 1))
                for i in range(self.Dimension - 1, self.current_bond, -1):
                    a = np.einsum('ijk,ljm,lk->mi', self.center[j][i], self.center[j][i].T, a)
                if (a-np.identity(a.shape[0])>0.01).any():
                    print('right_a:', a)
        else:
            for j in range(self.DataClass):
                a = np.ones((1, 1))
                for i in range(self.Dimension - 1, self.current_bond+1, -1):
                    a = np.einsum('ijk,ljm,lk->mi', self.center[j][i], self.center[j][i].T, a)
                if (a-np.identity(a.shape[0])>0.01).any():
                    print('right_a:', a)
    def judge_norm(self):
        #print('shape:',self.center[0][self.Dimension-1].shape)
        for j in range(self.DataClass):
            a = np.ones((1, 1))
            for i in range(self.Dimension):
                a = np.einsum('mi,ijk,njm->nk', a, self.center[j][i], self.center[j][i].T)
            if a-1>=1e-5:
                print('norm:',a)

    def judge_cumnlants(self):
        #k=self.current_bond
        for i in range(self.Number):
            for j in range(self.DataClass):
                a=np.ones((1,1))
                for k in range(self.current_bond):
                    if k!=0:
                        a=np.einsum('ij,jkl,k->il',a,self.center[j][k-1],self.tensor[i][k-1])

                        if (a-self.cumulants[i][j][k]>0.001).any():
                            print('i,j,k',i,j,k)


        for i in range(self.Number):
            for j in range(self.DataClass):
                a=np.ones((1,1))
                for k in range(self.Dimension-1,self.current_bond,-1):
                    if k!=self.Dimension-1:
                        a=np.einsum('ijk,j,kl->il',self.center[j][k+1],self.tensor[i][k+1],a)
                        if (a-self.cumulants[i][j][k]>0.001).any():
                            print('i,j,k',i,j,k)

    def judge_edge(self):
        for i in range(self.Number):
            for j in range(self.DataClass):

                if self.cumulants[i][j][0]!=1:
                    print('i,j,k', i, j,' left')

        for i in range(self.Number):
            for j in range(self.DataClass):

                if self.cumulants[i][j][self.Dimension-1]!=1:
                    print('i,j,k', i, j,' right')
    def init_cumnlants(self):


        for i in range(self.Number):
            self.cumulants[i]=[]
            InnerPro=np.zeros(self.DataClass)
            for j in range(self.DataClass):
                self.cumulants[i].append([])
                leftpart=rightpart=[np.ones((1,1))]
                #self.cumulants[i][j]=[]

                for k in range (self.Dimension-2,0,-1):
                    rightpart=[einsum('j,ijk,kl->il',self.tensor[i][k+1],self.center[j][k+1],rightpart[0])]+rightpart

                self.cumulants[i][j]=leftpart+rightpart

                InnerPro[j] = abs(einsum('ijk,j,klm,l,mn->in',self.center[j][0],self.tensor[i][0],self.center[j][1]
                                         ,self.tensor[i][1],rightpart[0]))


            self.NumberClass[i]=np.argmax(InnerPro)


    def merge_bond(self):
        k = self.current_bond
        for i in  range(self.DataClass):
            self.merged_matrix[i] = np.einsum('ijk,klm->ijlm', self.center[i][k],
                                       self.center[i][(k + 1) % self.Dimension], order='C')

    def normalize(self):
        self.merged_matrix = [self.merged_matrix[i]/norm(self.merged_matrix[i]) for i in range(self.DataClass)]

    def rebuild_bond(self,Dir):#还没计算cut的情况
        k = self.current_bond
        kp1 = (k + 1) % self.Dimension
        #print('merge:', k,kp1 )
        for i in range(self.DataClass):
            assert self.merged_matrix[i] is not None


            U, s, V = svd(self.merged_matrix[i].reshape((self.bond_dimension[i][
                                                          (k - 1) % self.Dimension] * 2, 2 * self.bond_dimension[i][kp1])))

            if s[0] <= 0.:
                print(
                    'Error: At bond %d Merged_mat happens to be all-zero.\nPlease tune learning rate.' % self.current_bond)
                raise FloatingPointError('Merged_mat trained to all-zero')

            #s_eff = s[s > self.cutoff*s[0]]
            s_eff=s


            l=min(len(s_eff),self.maxbond)

            s = np.diag(s[:l])
            U = U[:, :l]
            V = V[:l, :]
            self.bond_dimension[i][k] = l#Update k-th bond dimension

            if Dir == 'r':
                V = s@ V
                V /= norm(V)
                print('update:',k)
            else:
                U = U@ s
                U /= norm(U)
                print('Update:',kp1)

            self.center[i][k] = U.reshape((self.bond_dimension[i][(k - 1) % self.Dimension],
                                          2, l))
            self.center[i][kp1] = V.reshape((l, 2, self.bond_dimension[i][kp1]))

        self.current_bond += 1 if Dir == 'r' else -1
            #self.merged_matrix = None


    def Update_cumnlants(self,Dir):
        k=self.current_bond
        #print('cu_bond:',k)
        if Dir == 'r':


            for i in range(self.Number):
                #InnerPro=np.zeros(self.DataClass)
                for j in  range(self.DataClass):
                    self.cumulants[i][j][k]=einsum('li,ijk,j->lk',self.cumulants[i][j][k-1],self.center[j][k-1],self.tensor[i][k-1])

        elif Dir =='l':


            for i in range(self.Number):
                for j in range(self.DataClass):
                    self.cumulants[i][j][k+1]=einsum('j,ijk,kl->il',self.tensor[i][k+2],self.center[j][k+2],self.cumulants[i][j][k+2])

        else:
            print('not correct direction')

    def Update_Class(self):
        k=self.current_bond
        kp1 = (k + 1) % self.Dimension
        for i in range(self.Number):
            InnerPro=np.zeros(self.DataClass)
            for j in range(self.DataClass):

                InnerPro[j]=abs(np.einsum('ij,jklm,mn,k,l->in',self.cumulants[i][j][k],self.merged_matrix[j],
                                          self.cumulants[i][j][kp1],self.tensor[i][k],self.tensor[i][kp1]))
            self.NumberClass[i]=np.argmax(InnerPro)

    def bond_train(self,Dir):
        Number_in_class = np.zeros(self.DataClass)
        RightValue = {}
        k = self.current_bond
        kp1 = (k + 1) % self.Dimension
        km1 = (k - 1) % self.Dimension
        for j in range(self.DataClass):
            dl = self.bond_dimension[j][km1]
            dr = self.bond_dimension[j][kp1]

            RightValue[j] = np.zeros((dl, 2, 2, dr))

        for i in range(self.Number):
            # print('shape2:',RightValue[self.NumberClass[i]].shape)
            RightValue[self.NumberClass[i]] += np.einsum('ij,jkl,lmn,np->pmki',
                                                         self.cumulants[i][self.NumberClass[i]][kp1],
                                                         self.tensor[i][kp1].reshape(1, 2, 1),
                                                         self.tensor[i][k].reshape(1, 2, 1),
                                                         self.cumulants[i][self.NumberClass[i]][k])

            Number_in_class[self.NumberClass[i]] += 1

        for i in range(self.DataClass):
            if Number_in_class[i] != 0:
                self.merged_matrix[i] = RightValue[i]
        '''
        for i in range(self.Number):
            for j in range(self.DataClass):
                    if self.cumulants[i][j][k].shape[1]!=self.bond_dimension[j][km1]:
                        print(self.cumulants[i][j][k].shape[1])
                        print(self.bond_dimension[j][dl])
        '''
        self.normalize()
        self.rebuild_bond(Dir)

        self.Update_cumnlants(Dir)

        #self.judge_cumnlants()
        #self.judge_edge()

        self.merge_bond()

        self.Update_Class()

    def UpdateCenter(self):
        self.right_cano()
        self.init_cumnlants()

        for loop in range(5):

            for bond in range(0,self.Dimension-2):
                Dir ='r'
                self.bond_train(Dir)
                print('dir:',Dir)
                print('curr_bond:',self.current_bond)
                self.judge_left_cano(Dir)
                self.judge_right_cano(Dir)
                self.judge_norm()

            for bond in range(self.Dimension-2,0,-1):
                Dir = 'l'
                self.bond_train(Dir)
                print('dir:',Dir)
                print('curr_bond:',self.current_bond)
                self.judge_right_cano(Dir)
                self.judge_left_cano(Dir)
                self.judge_norm()



    def Accuracy(self):
        self.UpdateCenter()
        a=b=0
        for i in range(self.Number):
            if self.NumberClass[i]==0:
                a+=1
            elif self.NumberClass[i]==1:
                b+=1
            else:
                print('error')
        if a>b:
            class0=0
            class1 =1
        else:
            class0=1
            class1 = 0
        print(class0)

        #print(self.NumberClass)
        df = pd.read_csv(r'C:\Users\dell\Documents\GitHub\Dataset\breast-w_csv.csv')
        count=0
        for i in range(self.Number):
            if df.values[i][9]=='benign' and self.NumberClass[i]==class0:
                count+=1
            elif df.values[i][9] =='malignant' and self.NumberClass[i]==class1:
                count+=1
        print(count/self.Number)

