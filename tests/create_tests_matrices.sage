# We are creating the file "auto_tests_matrices.h".
# Heavily using Sage.
# https://sagemath.org


skel_def=u"""
SNmatrix<TYPE,SIZE> testsMatrix_NAME_A()
{
    SNmatrix<TYPE,SIZE> S(0);

    POPULATE_A

    return S;
}

SNupperTriangularMatrix<TYPE,SIZE> testsMatrix_NAME_A_U()
{
    SNupperTriangularMatrix<TYPE,SIZE> S(0);

    POPULATE_U

    return S;
}
"""

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
