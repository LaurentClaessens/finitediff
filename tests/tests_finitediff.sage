# This small Sage script creates the matrix A
# and computes its 'U' part.
# This is dedicated to the test
# test_finitediff_U


# https://sagemath.org


size=100
h=3

def c(k):
    return k

A=matrix(QQ,size)

for i in range(1,size-1):
    A[i,i]=2/(h**2)-c(i)
    A[i-1,i]=-1/h**2
    A[i+1,i]=-1/h**2

A[0,0]=2/h**2-c(0)
A[0,1]=-1/h**2
A[size-2,size-1]=-1/h**2
A[size-1,size-1]=2/h**2-c(size)

U=A.LU()[2]
for i in range(0,size):
    for j in range(0,size):
        if U[i,j] != 0:
            print("ans_U.at({},{})={};".format(i,j,numerical_approx(U[i,j])))
