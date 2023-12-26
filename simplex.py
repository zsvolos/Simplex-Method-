import numpy as np 
from pprint import pprint

def PhaseI(A: np.array, b: np.array, c: np.array, B: np.array, B_inv: np.array, bvar: np.array, nbvar: np.array):
    m = np.shape(A)[0]
    n = np.shape(A)[1]
    itr = 0
    xbvar = np.dot(B_inv, b)
    run = True
    
    while run:    

        xbvar = np.dot(B_inv, b)
        cb = np.zeros((1,m))
        for i in range(0,m):
            cb[0,i] = c[bvar[i]]
             #cb is a flat vector?!

        p = np.dot(cb,B_inv)
        cjbars = np.array([])
        newvar = -1
        newvarloc = -1
        temp = 0
        for i in range(0,len(nbvar)):
            j = nbvar[i]
            if c[j] - np.dot(p, A[:,j]) < temp: 
                newvar = j 
                newvarloc = i
                temp = c[j] - np.dot(p, A[:,j])
        if temp >= 0:
            print("all Cj are greater than zero")
            temp = n
            #check if an artificial var remains in basis
            for i in bvar:
                if n - (i+1) < m: 
                    return((-1, A , b, B_inv, bvar, nbvar))
            A = np.delete(A,-m, 1)
            return((0, A , b, B_inv, bvar, nbvar))  # no artificial vars in basis. 
            
        u = np.dot(B_inv, A[:,newvar])
        temp = 0
        for i in u:
            if i > temp:
                temp = i
        if temp == 0:
            print("not all u are positive")
            return((1, A , b, B_inv, bvar, nbvar))
            break


        theta = 1000000
        index = -1
        for i in range(0,m):
            if u[i] > 0:
                 if xbvar[i]/u[i] < theta:
                    theta = xbvar[i]/u[i]
                    index = i

        for i in range(0, len(xbvar)):
            if i == index:
                xbvar[i] = theta
            else:
                xbvar[i] += -theta*u[i]
        
        B[:,index] = A[:,newvar]
        
        temp = bvar[index] 
        bvar[index] = newvar
        nbvar = np.delete(nbvar, newvarloc)   # replaces nbvar[newvarloc] = temp

        B_inv = np.linalg.inv(B)
        itr += 1
        if itr > 1000:
            run = False

        itr += 1
    ###END OF WHILE####




def Simplex(A: np.array, b: np.array, c: np.array, B: np.array, B_inv: np.array, bvar: np.array, nbvar: np.array, blands = False):
    m = np.shape(A)[0]
    n = np.shape(A)[1]
    itr = 0
    xbvar = np.dot(B_inv, b)
    run = True
    
    while run:    

        xbvar = np.dot(B_inv, b)
        cb = np.zeros((1,m))
        for i in range(0,m):
            cb[0,i] = c[bvar[i]]
             #cb is a flat vector?!

        p = np.dot(cb,B_inv)
        cjbars = np.array([])
        newvar = -1
        newvarloc = -1
        temp = 0
        for i in range(0,len(nbvar)):
            j = nbvar[i]
            if c[j] - np.dot(p, A[:,j]) < temp: 
                newvar = j 
                newvarloc = i
                temp = c[j] - np.dot(p, A[:,j])
        if temp >= 0:
            print("all Cj are greater than or equal to zero")
            return((1,np.dot(cb,xbvar)[0,0] , xbvar, p, B_inv, bvar, nbvar))
            run = False
            break

        u = np.dot(B_inv, A[:,newvar])
        temp = 0
        for i in u:
            if i > temp:
                temp = i
        if temp == 0:
            print("not all u are positive")
            return((2, '-inf', xbvar, p, B_inv, bvar, nbvar))
            break


        theta = 1000000
        index = -1
        for i in range(0,m):
            if u[i] > 0:
                 if xbvar[i]/u[i] < theta:
                    theta = xbvar[i]/u[i]
                    index = i


        #update xbvar
        for i in range(0, len(xbvar)):
            if i == index:
                xbvar[i] = theta
            else:
                xbvar[i] += -theta*u[i]
        B[:,index] = A[:,newvar]
        
        temp = bvar[index] 
        bvar[index] = newvar
        nbvar[newvarloc] = temp

        B_inv = np.linalg.inv(B)
        itr += 1
        if itr > 1000:
            run = False
    ###END OF WHILE####








###PROBLEM ONE###

A = np.array([[-3, 2, 1,0],[-2,1,0,1]])
b = np.array([[30],[12]])
c = np.array([-5,-7,0,0])
B = np.array([[1,0],[0,1]])
B_inv = np.linalg.inv(B)
bvar = np.array([2,3])
nbvar = np.array([0,1])
print('PROBLEM ONE')
pprint(Simplex(A, b, c, B, B_inv, bvar, nbvar))


print('PROBLEM TWO')
###PROBLEM TWO###
A = np.array([[1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
              [4,2,2,2,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
              [0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
              [1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
              [0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
              [1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
              [0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0],
              [0,1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0],
              [0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0],
              [0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1]
], dtype= float)

b = np.transpose(np.array([[7,8,3,1.8,.3,3.8,3.2,.5,0.5,0.4]]))
c = np.array([-60,-40,-30,-30,-15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
cphase1 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1])
B = np.array([[1,0,0,0,0,0,0,0,0,0],
              [0,1,0,0,0,0,0,0,0,0],
              [0,0,1,0,0,0,0,0,0,0],
              [0,0,0,1,0,0,0,0,0,0],
              [0,0,0,0,1,0,0,0,0,0],
              [0,0,0,0,0,1,0,0,0,0],
              [0,0,0,0,0,0,1,0,0,0],
              [0,0,0,0,0,0,0,1,0,0],
              [0,0,0,0,0,0,0,0,1,0],
              [0,0,0,0,0,0,0,0,0,1]
], dtype= float)
B_inv = np.linalg.inv(B)
bvar = np.array([i for i in range(15,25)])
nbvar = np.array([i for i in range(0,15)])

print("PROB. 2 PHASE I")
p1 = PhaseI(A, b, cphase1, B, B_inv, bvar, nbvar)
print('PROB. 2 PHASE II ')
p2 = Simplex(p1[1],p1[2],c, np.linalg.inv(p1[3]), p1[3],p1[4],p1[5])
pprint(p2)





print('PROBLEM THREE')
###PROBLEM THREE###
A = np.array([[1/4,-8,-1,9,1,0,0],[1/2,-12,-1/2,3,0,1,0],[0.0,0,1,0,0,0,1]])
b = np.array([[0],[0],[1.0]])
c = np.array([-3/4, 20, -1/2, 6, 0, 0 ,0])
B = np.array([[1.0,0.0,0.0],[0,1,0],[0,0,1]])
B_inv = np.linalg.inv(B)
bvar = np.array([4,5,6])
nbvar = np.array([0,1,2,3])
pprint(Simplex(A, b, c, B, B_inv, bvar, nbvar))
print("cycled and ended at 1000 iterations")
#without blands rule, it cycles endlessly. here I have capped the iterations at 1000







