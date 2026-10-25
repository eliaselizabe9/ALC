#labo 6
import numpy as np
import random

def arroba(A, B):
    # Multiplicación matricial sin @
    A = np.asarray(A)
    B = np.asarray(B)
    
    if A.ndim == 2 and B.ndim == 2:
        res = np.zeros((A.shape[0], B.shape[1]))
        for i in range(A.shape[0]):
            for j in range(B.shape[1]):
                res[i, j] = np.sum(A[i, :] * B[:, j])
        return res
    elif A.ndim == 2 and B.ndim == 1:
        res = np.zeros(A.shape[0])
        for i in range(A.shape[0]):
            res[i] = np.sum(A[i, :] * B)
        return res
    elif A.ndim == 1 and B.ndim == 1:
        return np.sum(A * B)


def f_A(A, v):
    """Aplica f_A(v) = Av / ||Av||2."""
    w_prime = arroba(A, v)
    norma = np.linalg.norm(w_prime)
    if norma > 0:
        return w_prime / norma
    return np.zeros_like(v)


def metpot2k(A, tol=1e-15, K=1000):
    """
    Método de la potencia aplicando f_A dos veces por iteración (Algoritmo 1).
    Retorna el autovector v, el autovalor lambda y la cantidad de iteraciones k.
    """
    n = A.shape[0]
    v = np.random.randn(n)
    v = v / np.linalg.norm(v)
    
    v_tilde = f_A(A, f_A(A, v))
    e = arroba(v_tilde, v)
    k = 0
    
    while np.abs(e - 1.0) > tol and k < K:
        v = v_tilde
        v_tilde = f_A(A, f_A(A, v))
        e = arroba(v_tilde, v)
        k += 1
        
    # Estimación del autovalor lambda = v_tilde^T * A * v_tilde
    Av = arroba(A, v_tilde)
    lam = arroba(v_tilde, Av)
    
    return v_tilde, lam, k


def diagRH(A, tol=1e-15, K=1000):
    """
    Diagonalización de matriz simétrica A mediante Método de la Potencia y 
    Reflectores de Householder (Algoritmo 2).
    
    Retorna:
        - S: Matriz de autovectores.
        - D: Matriz de autovalores en la diagonal.
        Tales que A = S @ D @ S.T.
        
    Si la matriz A no es simétrica, retorna None.
    """
    # 1. Verificación de matriz cuadrada y simétrica
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None
    if not np.allclose(A, A.T, atol=1e-8):
        return None
        
    n = A.shape[0]
    
    # Caso base n = 1
    if n == 1:
        return np.array([[1.0]]), np.array([[A[0, 0]]])
        
    # 2. Calcular el autovector v1 y autovalor dominante lambda1 con metpot2k
    v1, lambda1, _ = metpot2k(A, tol, K)
    v1 = v1 / np.linalg.norm(v1)
    
    # 3. Construir el reflector de Householder Hv1 = I - 2*(e1 - v1)(e1 - v1)^T / ||e1 - v1||^2
    e1 = np.zeros(n)
    e1[0] = 1.0
    u = e1 - v1
    norma_u = np.linalg.norm(u)
    
    if norma_u < 1e-12:
        Hv1 = np.eye(n)
    else:
        Hv1 = np.eye(n) - 2.0 * np.outer(u, u) / (norma_u**2)
        
    # Caso base n = 2
    if n == 2:
        S = Hv1
        HA = arroba(Hv1, A)
        D = arroba(HA, Hv1.T)
        return S, D
    else:
        # 4. Transformar B = Hv1 * A * Hv1^T
        HA = arroba(Hv1, A)
        B = arroba(HA, Hv1.T)
        
        # 5. Paso recursivo sobre el bloque inferior A_bar (submatriz (n-1)x(n-1))
        A_bar = B[1:, 1:]
        S_bar, D_bar = diagRH(A_bar, tol, K)
        
        # 6. Reconstruir D y S
        D = np.zeros((n, n))
        D[0, 0] = lambda1
        D[1:, 1:] = D_bar
        
        M = np.zeros((n, n))
        M[0, 0] = 1.0
        M[1:, 1:] = S_bar
        
        S = arroba(Hv1, M)
        return S, D



#labo 7
import numpy as np
import random 

def transiciones_al_azar_continuas(n):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    Retorna matriz T de n x n normalizada por columnas,
    y con entradas al azar en el intervalo [0,1]
    """
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            res[i][j] = random.uniform(0.0, 1.0)

    return res / np.sum(res, axis=0)

'''
otra forma sin usar axis:
def transiciones_al_azar_continuas(n):
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            res[i][j] = random.uniform(0.0, 1.0)

    for j in range(n):
        norma = 0
        for i in range(n):
            norma += res[i][j]
        for i in range(n):
            res[i][j] = res[i][j] / norma
    return res
'''

def transiciones_al_azar_uniformes(n, thres):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    thres probabilidad de que una entrada sea distinta de cero.
    Retorna matriz T de n x n normalizada por columnas. 
    El elemento i,j es distinto de cero si el número generado al azar para i,j es menor o igual a thres. 
    Todos los elementos de la columna $j$ son iguales 
    (a 1 sobre el número de elementos distintos de cero en la columna).
    """
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            res[i][j] = random.uniform(0.0, 1.0)

    for i in range(n):
        for j in range(n):
            if res[i][j] <= thres:
                res[i][j] = 1
            else: res[i][j] = 0

    for j in range(n):
        suma_col = 0
        for i in range(n):
            suma_col += res[i][j]
        if suma_col == 0:
            indice = random.randint(0, n-1)
            res[indice][j] = 1.0    

    return res / np.sum(res, axis=0)

def nucleo(A, tol=1e-15):
    """
    A una matriz de m x n
    tol la tolerancia para asumir que un vector esta en el nucleo.
    Calcula el nucleo de la matriz A diagonalizando la matriz traspuesta(A) * A (* la multiplicacion matricial), usando el medodo diagRH. El nucleo corresponde a los autovectores de autovalor con modulo <= tol.
    Retorna los autovectores en cuestion, como una matriz de n x k, con k el numero de autovectores en el nucleo.
    """
    m, n = A.shape
    B = arroba(A.T, A) # B de dim nxn
    S,D = diagRH(B, tol)

    indices = []
    k = 0
    for i in range(n):
        if np.abs(D[i][i]) <= tol:
            indices.append(i)
            k += 1

    if k == 0:
        return np.zeros((0, 0))

    res = np.zeros((n, k))
    for j in range(k):
        l = indices[j]
        for i in range(n):
            res[i][j] = S[i][l]
    return res

def crea_rala(listado, m_filas, n_columnas, tol=1e-15):
    """
    Recibe una lista listado, con tres elementos: lista con indices i, lista con indices j, y lista con valores A_ij de la matriz A. 
    Tambien las dimensiones de la matriz a traves de m_filas y n_columnas. Los elementos menores a tol se descartan.
    Idealmente, el listado debe incluir unicamente posiciones correspondientes a valores distintos de cero. Retorna una lista con:
    - Diccionario {(i,j):A_ij} que representa los elementos no nulos de la matriz A. Los elementos con modulo menor a tol deben descartarse por default. 
    - Tupla (m_filas,n_columnas) que permita conocer las dimensiones de la matriz.
    """
    dim = (m_filas, n_columnas)
    dic = dict()

    if len(listado) == 0 or len(listado[0]) == 0: 
        return [dic, dim] 
    
    indices_i = listado[0]
    indices_j = listado[1]
    valores_Aij = listado[2]
    for k in range(len(indices_i)):
        if abs(valores_Aij[k]) > tol:
            dic[indices_i[k], indices_j[k]] = valores_Aij[k]
    return [dic, dim]

def multiplica_rala_vector(A,v):
    """
    Recibe una matriz rala creada con crea_rala y un vector v. 
    Retorna un vector w resultado de multiplicar A con v
    """
    m, n = A[1]
    v = np.asarray(v)

    if n != len(v):
        return

    res = np.zeros(m)
    for i,j in A[0]:
        res[i] += A[0][i,j] * v[j]

    return res
