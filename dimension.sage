from symbolictetra.py import complement, doublecomplement

def WD_matrix(lengths):
    """
    Parameters
    List of Edgelengths [12,13,14,23,24,34]
    
    Returns:
    The D matrix!
    """
    D=Matrix([[0, lengths[0]**2, lengths[1]**2, lengths[2]**2, 1], 
             [lengths[0]**2, 0, lengths[3]**2, lengths[4]**2, 1],
             [lengths[1]**2, lengths[3]**2, 0, lengths[5]**2, 1],
             [lengths[2]**2, lengths[4]**2,lengths[5]**2,0,1],
             [1,1,1,1,0]])
    return D


def lexi(arr):
    return [arr[0], arr[2], arr[4], arr[5], arr[3], arr[1]]


def D_3(i,j,k, lengths):
    """
    Parameters
    lengths:List of Edgelengths [12,13,14,23,24,34]; one or more of which is a symbol (variable).
    Distinct i,j,k in {1,2,3,4}
    Returns:
    D_ijk, as defined in WD's paper
    """
    l=complement(i, j, k)
    M=WD_matrix(lengths)
    M=M.delete_columns([l-1])
    M=M.delete_rows([l-1])
    return M.det()

def D_2(i,j, lengths):
    """
    Parameters
    i,j: distinct integers in {1,2,3,4}
    lengths: List of Edgelengths [12,13,14,23,24,34]; one or more of which is a symbol (variable).
    Returns:
    D_ij, as defined in WD's paper
    """
    assert i!=j and i in [1,2,3,4] and j in [1,2,3,4]
    k,l=doublecomplement(i, j)
    M=WD_matrix(lengths)
    M=M.delete_columns([k-1])
    M=M.delete_rows([l-1])
    return (-1)**(k+l)*M.det()

def lengths_to_dihedral_exponentials(lengths):
    """
    Parameters: 
    lengths: List of Edgelengths [12,13,14,23,24,34]
    """
    angles={}
    for i in range(1,5):
        for j in range(1,5):
            if i<j:
                k,l=doublecomplement(i, j)
                real_part=D_2(i,j,lengths)/sqrt( D_3(i,j,k, lengths)* D_3(i,j,l, lengths))
                imaginary_part=sqrt(1-real_part**2)
                angles[(i,j)]=real_part**2+2*real_part*imaginary_part*I-imaginary_part**2
    return angles

def lengths_to_dihedral_exponential_single(lengths):
    """
    Parameters: 
    lengths: List of Edgelengths [12,13,14,23,24,34]
    """
    angles={}
    for i in range(1,5):
        for j in range(1,5):
            if i<j:
                k,l=doublecomplement(i, j)
                real_part=D_2(i,j,lengths)/sqrt( D_3(i,j,k, lengths)* D_3(i,j,l, lengths))
                imaginary_part=sqrt(1-real_part**2)
                angles[(i,j)]=real_part+imaginary_part*I
    return angles

def lengths_to_dihedrals(lengths):
    """
    Parameters: 
    lengths: List of Edgelengths [12,13,14,23,24,34]
    """
    angles={}
    for i in range(1,5):
        for j in range(1,5):
            if i<j:
                k,l=doublecomplement(i, j)
                v=D_2(i,j,lengths)/sqrt( D_3(i,j,k, lengths)* D_3(i,j,l, lengths))
                angles[(i,j)]=arccos(v,hold=True)
    return angles

def dimensions(v):
    exponentials=lengths_to_dihedral_exponentials(v)
    D=WD_matrix(v).det()
#     print('the value of D is', D)
    K = QuadraticField(-2*D)
    L = K.ring_of_integers()
    A=[]
    
    primes=set()
    for edge in exponentials:
        for p,e in L.fractional_ideal(exponentials[edge]).factor():
            primes.add(p)
            
            
            
    for p in primes:
        p_valuations=[]
        for edge in [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]:
            valuation=0
#             print(edge, L.fractional_ideal(exponentials[edge]).factor())
            for q,f in L.fractional_ideal(exponentials[edge]).factor():
                if q==p:
                    valuation=f
            p_valuations.append(valuation)
        A.append(p_valuations)
    val_matrix=matrix(A)
    return val_matrix.rank()

def check_dehn_invariant(vec):
    D=WD_matrix(vec).det()
    
    if not is_square(2*D) and not is_square(6*D):
        return False       
    
    #input is in kkpr format
    value = product([lengths_to_dihedral_exponentials(vec)[edge] for edge in [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]])
    if bool(value**4==1) or bool(value**3==1):
        return True
    return False


def kimi_family(t):
    assert t>4
    return [t+1, t+1, t, 6, t-1, t-1]

def hill_family_1(a,b):
    x=3*b**2-a**2
    y=6*a*b
    return [x,x,y,x,a**2+3*b**2,a**2+3*b**2]



