import numpy as np
import torch
from torch import einsum
from scipy.linalg import norm
import math
import tensorflow.keras as keras
import pandas as pd

class Update:
    def __init__(self, tensor, center, Number, DataClass, Dimension, bond_dimension, maxbond, cutoff=10**-9):

        self.tensor = tensor  # 数据点的mps
        self.center = center  # 中心点的mps
        self.Number = Number  # 一共有Number个数
        self.DataClass = DataClass  # 将数据分为DataClass个类别
        self.Dimension = Dimension  # 数据一共有多少维
        self.current_bond = 0  # 当前更新的bond编号
        self.cumulants = [torch.ones(self.Number,1).double() for i in range(self.DataClass)]   #存储
        #self.NumberClass = {}  # 数据类别
        self.bond_dimension = []
        for i in range(self.DataClass):
            self.bond_dimension.append(list(bond_dimension))  # 中心点之间bond维数

        self.merged_matrix = [None for i in range(self.DataClass)]  # 合并起来的4阶张量
        self.maxbond = maxbond  # 最大的bond维数
        self.cutoff = cutoff


    def right_cano(self):

        for i in range(self.DataClass):
            for j in range(self.Dimension - 1, -1, -1):
                dl = self.center[i][j].shape[0]
                dr = self.center[i][j].shape[2]
                U, s, V = torch.svd(self.center[i][j].reshape(dl, dr * 2).T)
                myd = min(dl, dr * 2)
                Q = U[:, :myd]
                # R=np.diag(s)@(V[:,:myd].T)
                self.center[i][j] = Q.T.reshape(-1, 2, dr)

    def init_cumnlants(self):

        #print('center:',self.center[i][n])

        for i in range(self.DataClass):
            rightpart = [torch.ones((1, self.Number)).double()]
            for n in range(self.Dimension - 1, 1, -1):
                #print('shape:',self.center[i][n].shape,rightpart[0].shape)
                rightpart = [einsum('ij,ljr,ri -> li',
                                    self.tensor[:, n, :],
                                    self.center[i][n],
                                    rightpart[0])] + rightpart

            self.cumulants[i] = [self.cumulants[i]] + rightpart


        Storage = np.zeros((self.Number,self.DataClass))
        for i in range(self.DataClass):
            Storage[:, i] = einsum('jkl,lmn,ik,im,ni->ji',
                                   self.center[i][0],self.center[i][1],
                                   self.tensor[:,0,:],
                                   self.tensor[:,1,:],
                                   self.cumulants[i][1]).squeeze()



        Storage = np.abs(Storage)
        self.NumberClass = [np.argmax(c) for c in Storage]
        #print('class:',self.NumberClass)

    def merge_bond(self):
        k = self.current_bond
        for i in range(self.DataClass):

            self.merged_matrix[i] = einsum('ijk,klm->ijlm', self.center[i][k],self.center[i][(k + 1) % self.Dimension])

    def normalize(self):
        for i in range(self.DataClass):
            self.merged_matrix[i] /= norm(self.merged_matrix[i])

    def rebuild_bond(self, Dir):
        k = self.current_bond
        kp1 = (k + 1) % self.Dimension
        # print('merge:', k,kp1 )
        for i in range(self.DataClass):
            assert self.merged_matrix[i] is not None

            U, s, V = torch.svd(self.merged_matrix[i].reshape((self.bond_dimension[i][
                                                             (k - 1) % self.Dimension] * 2,
                                                         2 * self.bond_dimension[i][kp1])))


            s_eff = s[s > self.cutoff * s[0]]

            if len(s_eff) < 2:
                s_eff = s[:2]

            l = min(len(s_eff), self.maxbond)

            s = torch.diag(s[:l])
            U = U[:, :l]
            V = V.T[:l, :]#numpy里u@s@v=matrix,torch里u@s@v.T=matrix

            self.bond_dimension[i][k] = l  # Update k-th bond dimension

            if Dir == 'r':
                V = s @ V
                V /= norm(V)
                # print('update:',k)
            else:
                U = U @ s
                U /= norm(U)
                # print('Update:',kp1)
            #print('before:',self.merged_matrix[0])
            self.center[i][k] = torch.reshape(U,(self.bond_dimension[i][(k - 1) % self.Dimension], 2, l))

            self.center[i][kp1] =torch.reshape(V,(l, 2, self.bond_dimension[i][kp1]))

            # print('i',i,'bond_d:', self.bond_dimension)
        # print('bond_d:', self.bond_dimension[0][k])
        print('update ', k, '-th site...')
        self.current_bond += 1 if Dir == 'r' else -1
        # self.merged_matrix = None

        # print('rebulid:',self.bond_dimension)

    def Update_cumnlants(self, Dir):
        k = self.current_bond
        # print('cu_bond:',k)

        if Dir == 'r':
            for i in range(self.DataClass):
                self.cumulants[i][k] = einsum('ij,jlk,il->ik',self.cumulants[i][k-1],
                                              self.center[i][k-1],
                                              self.tensor[:,k-1,:])


        elif Dir == 'l':
            for i in range(self.DataClass):
                self.cumulants[i][k+1] = einsum('jl,ilk,kj->ij',self.tensor[:,k+2,:],self.center[i][k+2],self.cumulants[i][k+2])

    def judge(self):
        k=self.current_bond
        kp1=self.current_bond+1
        #是否左规范
        for i in range(self.DataClass):
            leftpart = torch.ones(1,1,1,1).double()
            #print('center shape:',self.center[i][j].shape)
            for j in range(k):
                leftpart = einsum('ijkl,jmn,lma->inka',leftpart,self.center[i][j],self.center[i][j])
            leftpart = leftpart.squeeze()
            leftpart = leftpart.numpy()
            a = np.identity(self.center[i][j].shape[2])
            m = leftpart - a
            m=np.abs(m)
            '''
            flag =True
            for row in range(self.center[i][j].shape[2]):
                for column in range(self.center[i][j].shape[2]):
                    if m[row][column] > 10**-7:
                        flag = False
                        break
            
            if flag == True:
                print('left cano')
            else:
                print('not left cano')
            '''
            print('m1:',m)
        #是否右规范化
        for i in range(self.DataClass):
            rightpart = torch.ones(1,1,1,1).double()
            for j in range(self.Dimension-1,kp1,-1):
                #print('center shape:',self.center[i][j].shape,rightpart.shape)
                rightpart = einsum('jnla,imj,kml->inka',rightpart,self.center[i][j],self.center[i][j])
            rightpart = rightpart.squeeze()
            rightpart = rightpart.numpy()
            a = np.identity(self.center[i][j].shape[0])
            m = rightpart - a
            m=np.abs(m)
            '''
            flag = True
            for row in range(self.center[i][j].shape[0]):
                for column in range(self.center[i][j].shape[0]):
                    if m[row][column] > 10**-7:
                        flag = False
                        break
            
            if flag == True:
                print('right cano')
            else:
                print('not right cano')
            '''
            print('m2:', m)
        #是否内积为1
        for i in range(self.DataClass):
            leftpart = torch.ones(1,1,1,1).double()
            for j in range(self.Dimension):
                leftpart = einsum('ijkl,jmn,lma->inka',leftpart,self.center[i][j],self.center[i][j])
            print('Inner product:',leftpart.squeeze())
    def Update_Class(self):
        k = self.current_bond
        kp1 = (k + 1) % self.Dimension


        #Distance=torch.zeros((self.Number,self.DataClass))
        Storage = np.zeros((self.Number, self.DataClass))
        for i in range(self.DataClass):
            Storage[:, i] = einsum('ij,jklm,ik,il,mi->i', self.cumulants[i][k],
                                   self.merged_matrix[i],
                                   self.tensor[:, k, :],
                                   self.tensor[:, kp1, :],
                                   self.cumulants[i][kp1])

        #a=self.NumberClass
        Storage = np.abs(Storage)
        self.NumberClass = [np.argmax(c) for c in Storage]
        b = [np.max(c) for c in Storage]
        print('quantum cost:', sum(b))
        cost = 0
        for i in range(self.Number):
            cost += np.max(Storage[i,:])
        print('cost:',cost)

    def bond_train(self, Dir):
        Number_in_class = np.zeros(self.DataClass)
        RightValue = {}
        k = self.current_bond
        kp1 = (k + 1) % self.Dimension
        km1 = (k - 1) % self.Dimension
        for j in range(self.DataClass):
            dl = self.bond_dimension[j][km1]
            dr = self.bond_dimension[j][kp1]

            RightValue[j] = torch.zeros(dl, 2, 2, dr).double()

        for i in range(self.Number):

            RightValue[self.NumberClass[i]] += einsum('i,j,k,l->lkji',
                                                      self.cumulants[self.NumberClass[i]][kp1][:,i],
                                                      self.tensor[i][kp1],
                                                      self.tensor[i][k],
                                                      self.cumulants[self.NumberClass[i]][k][i,:])

            Number_in_class[self.NumberClass[i]] += 1

        for i in range(self.DataClass):
            if Number_in_class[i] != 0:
                self.merged_matrix[i] = RightValue[i]
            else:
                self.merged_matrix[i] = einsum('ijk,klm->ijlm', self.center[i][k], self.center[i][kp1])

        #print('merge matrix:',self.merged_matrix[0])

        self.normalize()

        self.rebuild_bond(Dir)

        self.Update_cumnlants(Dir)

        self.merge_bond()

        self.Update_Class()
        # print('bond_dim:',self.bond_dimension)

    def UpdateCenter(self):

        self.right_cano()

        self.init_cumnlants()


        for bond in range(0, self.Dimension - 2):
            Dir = 'r'
            self.bond_train(Dir)
            self.judge()


        for bond in range(self.Dimension - 2, 0, -1):
            Dir = 'l'
            self.bond_train(Dir)
            self.judge()
        
        print('bond dimension:', self.bond_dimension)
        (x_test, y_test) = keras.datasets.mnist.load_data()
        train_labels = x_test[1]  # mnist.test.labels

        counts = np.zeros((self.DataClass, 10))
        for i in range(self.Number):
            counts[self.NumberClass[i]][train_labels[i]] += 1

        # 将最频繁的标签分配给质心
        a = [np.max(c) for c in counts]
        print('训练集准确率:', sum(a) / self.Number)

    def Accuracy(self):  # 1.test-number需要修改
        self.UpdateCenter()
        print('bond dimension:', self.bond_dimension)

        (x_test, y_test) = keras.datasets.mnist.load_data()
        # 导入MNIST数据集
        TestNumber = 1000
        test_data = {}  # mnist.test.images
        for i in range(TestNumber):
            test_data[i] = y_test[0][i].reshape(28 * 28) / 255

        TestTensor = np.zeros((TestNumber,self.Dimension,2))
        for i in range(TestNumber):
            for j in range(self.Dimension):
                TestTensor[i][j] = [math.cos(test_data[i][j]*math.pi/2),math.sin(test_data[i][j]*math.pi/2)]

        TestTensor = torch.from_numpy(TestTensor)

        train_labels = x_test[1]  # mnist.test.labels
        test_labels = y_test[1]  # mnist.train.labels
        counts = np.zeros((self.DataClass, 10))
        for i in range(self.Number):
            counts[self.NumberClass[i]][train_labels[i]] += 1
        labels_map = [np.argmax(c) for c in counts]
        # 将最频繁的标签分配给质心
        a = [np.max(c) for c in counts]
        print('训练集准确率:', sum(a) / self.Number)


        Storage = np.zeros((TestNumber,self.DataClass))
        for i in range(self.DataClass):
            rightpart = torch.ones((1, TestNumber)).double()
            for n in range(self.Dimension - 1, -1, -1):
                rightpart = einsum('ij,ljk,ki->li',
                                         TestTensor[:, n, :],
                                         self.center[i][n],
                                         rightpart)
            Storage[:, i] = rightpart.squeeze()

        TestNumberClass = [np.argmax(c) for c in Storage]



        Correct = 0
        for i in range(TestNumber):
            if labels_map[TestNumberClass[i]] == test_labels[i]:
                Correct += 1

        print("测试集准确率为: %f " % ((Correct / TestNumber) * 100))





