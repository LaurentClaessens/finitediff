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

SNpermutation<TYPE,SIZE> testMatrix_NAME_A_P()
{
    SNpermutation<TYPE,SIZE> S(0);

    POPULATE_P

    return S;
}

SNlowerTriangular<TYPE,SIZE> testsMatrix_NAME_A_L()
{
    SNlowerTriangular<TYPE,SIZE> S(0);

    POPULATE_L

    return S;
}

SNupperTriangular<TYPE,SIZE> testsMatrix_NAME_A_U()
{
    SNupperTriangular<TYPE,SIZE> S(0);

    POPULATE_U

    return S;

}
"""

def populate(A):
    """
    - A is a matrix
    Return the C++ code that populates a SNmatrix with the
      non vanishing elements of A
    """
    text=[]
    size=A.dimensions()[0]
    for i in range(0,size):
        for j in range(0,size):
            if A[i,j] != 0:
                text.append("S.at({},{})={};".format(i,j,numerical_approx(A[i,j])))
    return "\n".join(text)


def matrix_to_cpp(A,ttype,name):
    # - A : a matrix
    # - ttype : the (C++) type of the matrix elements (typically : "double")
    #    (this is a string)
    # - name : (string), the name of the matrix
    # return the C++ code that creates the matrix and its P,L,U.

    plu=A.LU()
    P=plu[0]
    L=plu[1]
    U=plu[2]
    new_text=skel_def
    size=str(A.dimensions()[0])
    new_text=new_text.replace("SIZE",size)
    new_text=new_text.replace("TYPE",ttype)
    new_text=new_text.replace("NAME",name)
    new_text=new_text.replace("POPULATE_A",populate(A))
    new_text=new_text.replace("POPULATE_P",populate(P))
    new_text=new_text.replace("POPULATE_L",populate(L))
    new_text=new_text.replace("POPULATE_U",populate(U))
    return new_text


A=matrix( [ [1,2],[3,4] ]  )


print(matrix_to_cpp(A,"int","FOO"))
