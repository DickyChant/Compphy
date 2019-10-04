#!/usr/bin/env python
# coding: utf-8

# # Define Class: vector and square matrix

# In[1]:


class vector:
    N = 0
    Vec = []

    def __init__(self, vec=None):
        if vec == None or vec == []:
            self.N = 0
            self.Vec = []
        else:
            self.Vec = [i for i in vec]
            self.N = len(vec)

    def print_vec(self):
        print(self.Vec)
    def print_vec_procession(self,p = 4):
        print("(",end='')
        for i in range(self.N):
            if(i==0):
                print(" %.*f"%(p,self.Vec[i]),end='')
            else:
                print(", %.*f"%(p,self.Vec[i]),end='')
        print(" )")


# In[2]:


class square:
    N = 0
    Squ = []

    def __init__(self, n=0, squ=None):
        if squ == None or squ == []:
            self.N = 0
            self.Squ = []
        else:
            self.N = n
            self.Squ = [[j for j in i] for i in squ]

    def print_squ(self):
        for row in self.Squ:
            print(row)

    def transpose(self):
        temp = [[self.Squ[j][i] for j in range(self.N)] for i in range(self.N)]
        return square(self.N, temp)


# # Define function: matrix multiply matrix and matrix multiply vector

# In[3]:


def matmul(m1, m2):
    if m1.N != m2.N:
        print("Error")
        return -1
    else:
        tmp1 = []
        for i in range(m1.N):
            tmp2 = []
            for j in range(m1.N):
                tmp3 = 0
                for k in range(m1.N):
                    tmp3 += m1.Squ[i][k]*m2.Squ[k][j]
                tmp2.append(tmp3)
            tmp1.append(tmp2)
        return square(m1.N, tmp1)


# In[4]:


def mvmul(mat, vec):
    if mat.N != vec.N:
        print("Error")
        return -1
    else:
        tmp1 = []
        for i in range(mat.N):
            tmp2 = 0
            for j in range(mat.N):
                tmp2 += mat.Squ[i][j]*vec.Vec[j]
            tmp1.append(tmp2)
        return vector(tmp1)


# # Hilbert Matrix Generator

# In[5]:


def Hilbert_mat(n):
    mat = [[1/(i+j+1) for j in range(n)] for i in range(n)]
    return square(n, mat)


# # GEM code

# In[6]:


def GEM(mat, vec=None):
    if vec != None and mat.N != vec.N:
        print("Error")
        return -1
    if vec == None:
        vec = vector([1 for i in range(mat.N)])
    mat0 = square(mat.N, mat.Squ)
    vec0 = vector(vec.Vec)
    for i in range(mat0.N):

        maxi = i
        for j in range(mat0.N-i):
            if abs(mat0.Squ[i+j][i]) > abs(mat0.Squ[maxi][i]):
                maxi = j+i
        if mat0.Squ[maxi][i] == 0:
            return -1
        temp = mat0.Squ[maxi]
        mat0.Squ[maxi] = mat0.Squ[i]
        mat0.Squ[i] = temp
        temp = vec0.Vec[maxi]
        vec0.Vec[maxi] = vec0.Vec[i]
        vec0.Vec[i] = temp

        temp = -1.0/mat0.Squ[i][i]
        for j in range(mat.N-i-1):
            l = temp*mat0.Squ[i+j+1][i]
            for k in range(mat.N-i):
                mat0.Squ[i+j+1][i+k] += l*mat0.Squ[i][i+k]
            vec0.Vec[i+j+1] += l*vec0.Vec[i]

    return mat0, vec0


# # Square root: using Newton method

# In[7]:


def sqr_newton(s, n=20):
    x = s
    for i in range(n):
        x = x/2+s/2/x

    return x


# # Cholesky code

# In[8]:


def Cholesky(mat):
    temp = [[0 for i in range(mat.N)] for j in range(mat.N)]
    for i in range(mat.N):
        if i == 0:
            temp[0][0] = sqr_newton(mat.Squ[0][0])
            continue
        for j in range(i):
            temp[i][j] = mat.Squ[i][j]
            for k in range(j):
                temp[i][j] -= temp[i][k]*temp[j][k]
            temp[i][j] /= temp[j][j]
        temp[i][i] = mat.Squ[i][i]
        for j in range(i):
            temp[i][i] -= temp[i][j]**2
        temp[i][i] = sqr_newton(temp[i][i])

    return square(mat.N, temp)


# # upper&lower triangular matrix solver

# In[9]:


def upper_tri_solver(mat, vec=None):
    if mat.N != vec.N and vec != None:
        print("Error")
        return -1
    if vec == None:
        vec = vector([1 for i in range(mat.N)])

    temp_vec = []
    for i in range(mat.N):
        if mat.Squ[i][i] == 0:
            print("Cannot Solve!")
            return -1
        temp = vec.Vec[mat.N-1-i]
        for k in range(i):
            temp -= mat.Squ[mat.N-1-i][mat.N-1-k]*temp_vec[k]
        temp /= mat.Squ[mat.N-1-i][mat.N-1-i]
        temp_vec.append(temp)
    temp_vec.reverse()
    return vector(temp_vec)


# In[10]:


def lower_tri_solver(mat, vec=None):
    if mat.N != vec.N and vec != None:
        print("Error")
        return -1
    if vec == None:
        vec = vector([1 for i in range(mat.N)])

    temp_vec = []
    for i in range(mat.N):
        if mat.Squ[i][i] == 0:
            print("Cannot Solve!")
            return -1
        temp = vec.Vec[i]
        for k in range(i):
            temp -= mat.Squ[i][k]*temp_vec[k]
        temp /= mat.Squ[i][i]
        temp_vec.append(temp)
    return vector(temp_vec)


# # Results

# In[11]:


for i in range(10):
    print("------n={}------".format(i+1))

    n = i+1
    b = vector([1 for i in range(n)])
    H = Hilbert_mat(n)

    print("------GEM------")

    H_GEM, b_GEM = GEM(H, b)
    x_GEM = upper_tri_solver(H_GEM,b_GEM)
    x_GEM.print_vec_procession()
    
    print("------Cholesky------")

    U = Cholesky(H)
    U_dag = U.transpose()
    y = lower_tri_solver(U,b)
    x_Cho = upper_tri_solver(U_dag,y)
    x_Cho.print_vec_procession()

