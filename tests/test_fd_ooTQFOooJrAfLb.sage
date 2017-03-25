# This small Sage script creates the matrix A
# and computes its 'U' part.
# This is dedicated to the test
# test_finitediff_U


# https://sagemath.org

# We create the C++ file
# ooTQFOooJrAfLb.h
# which contains the code for the test matrix and its answer.

size=100
TYPE="double"
SIZE=str(size)
h=3

def c(k):
    return k

def populate(A):
    """
    - A is a matrix
    - we write the C++ code that populates a SNmatrix with the
      non vanishing elements of A
    """
    text=[]
    for i in range(0,size):
        for j in range(0,size):
            if A[i,j] != 0:
                text.append("S.at({},{})={};".format(i,j,numerical_approx(A[i,j])))
    return "\n".join(text)


skel=open("ooTQFOooJrAfLb.skel.h","r").read()
print(skel)

def fd_matrix(h,size):
    A=matrix(QQ,size)

    for i in range(1,size-1):
        A[i,i]=2/(h**2)-c(i)
        A[i-1,i]=-1/h**2
        A[i+1,i]=-1/h**2

    A[0,0]=2/h**2-c(0)
    A[0,1]=-1/h**2
    A[size-2,size-1]=-1/h**2
    A[size-1,size-1]=2/h**2-c(size)
    
    return A

A=fd_matrix(h,size)
U=A.LU()[2]

new_text=skel.replace("TYPE",TYPE).replace("SIZE",SIZE)
new_text=new_text.replace("POPULATE_A",populate(A))
new_text=new_text.replace("POPULATE_U",populate(U))

print(new_text)
open("ooTQFOooJrAfLb.h","w").write(new_text)
