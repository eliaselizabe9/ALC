import numpy as np

def f(A, v):
    """f_A(v) = Av / ||Av||"""
    w = A @ v
    return w / norma2(w)

def f2(A, v):
    """f_A^2(v) = f_A(f_A(v))"""
    return f(A, f(A, v))


def norma2(v): #calcula la norma 2 de un vector
    v = np.array(v)
    n = len(v)
    contador = 0

    for i in range(n):
        contador+= v[i]**2

    return contador**(1/2)

def metpot2k(A, tol=1e-15, K=1000):
    '''A una matriz de nxn.
    tol la tolerancia en la diferencia entre un paso y el siguiente de la estimacion del autovector.
    K el numero maximo de iteraciones a realizarse.
    Retorna vector v, autovalor lambda y numero de iteraciones realizadas k'''
    n = A.shape[0]
    v = np.random.randn(n) #busca un vector al azar de tamano n
    v = v / norma2(v)  #normaliza v 
    v_t = f2(A, v) # hace la mutiliplacacion matricial y normaliza
    e = v_t @ v
    k = 0
    while abs(e - 1) > tol and k < K:
        v = v_t
        v_t = f2(A, v)
        e = v_t @ v
        k += 1
    lam = v_t @ A @ v_t            # cociente de Rayleigh (||v_t|| = 1)
    return v_t, lam, k

def diagRH(A, tol=1e-15, K=1000):
    n = A.shape[0]
    if n == 1:
        return np.eye(1), A.copy()

    v1, l1, _ = metpot2k(A, tol, K)
    e1 = np.zeros(n)
    e1[0] = 1.0

    diff = e1 - v1
    norm_diff = np.sqrt(diff @ diff)
    if norm_diff > 1e-15:
        u = diff / norm_diff
        H_v1 = np.eye(n) - 2.0 * np.outer(u, u)
    else:
        H_v1 = np.eye(n)

    # Transformacion ortogonal por bloques
    B = np.matmul(np.matmul(H_v1, A), H_v1)

    A_tilde = B[1:, 1:]
    S_tilde, D_tilde = diagRH(A_tilde, tol, K)

    D = np.zeros((n, n))
    D[0, 0] = l1
    D[1:, 1:] = D_tilde

    S_aux = np.eye(n)
    S_aux[1:, 1:] = S_tilde

    S = np.matmul(H_v1, S_aux)
    return S, D




def _calcular_bases_svd(A, tol=1e-15):
    """
    Función interna que evalúa dimensiones, diagonaliza y calcula 
    los componentes base de la SVD (valores singulares, V y B).
    """
    m, n = A.shape
    transpuesta = False
    
    # Optimización por dimensiones (m < n vs m > n)
    if m < n:
        A = A.T
        transpuesta = True
        
    # Diagonalización de A^t * A
    AtA = np.dot(A.T, A) #producto entre A y A^t
    autovalores, autovectores = diagRH(AtA)  # obtengo autovalores y autovectores de A * At
    
    # Ordenar de mayor a menor
    idx_ordenados = np.argsort(autovalores)[::-1] 
    autovalores = autovalores[idx_ordenados]
    autovectores = autovectores[:, idx_ordenados]
    
    # Calcular valores singulares
    valores_singulares = np.sqrt(np.maximum(autovalores, 0))
    
    # Rango efectivo (valores mayores a la tolerancia)
    r = np.sum(valores_singulares > tol)
    
    # Calcular B = A * V para luego obtener U
    B = np.dot(A, autovectores)
    
    if transpuesta:
        # Si se traspuso, intercambiamos roles para devolver U y V correctamente ordenados
        return B, valores_singulares, autovectores, r, True
    else:
        return autovectores, valores_singulares, B, r, False


def svd_reducida(A, k="max", tol=1e-15):
    """
    Calcula la SVD reducida de la matriz A (de m x n).
    Parámetros:
    - A: matriz de m x n
    - k: número de valores singulares a retener ("max" o un entero)
    - tol: tolerancia para valores nulos
    Retorna:
    - hatU: matriz de m x k con los vectores singulares izquierdos
    - hatSig: vector de k valores singulares
    - hatV: matriz de n x k con los vectores singulares derechos
    """
    m, n = A.shape
    
    # Obtenemos los componentes base de la función auxiliar
    res1, valores_singulares, res3, r, transpuesta = _calcular_bases_svd(A, tol)
    
    # Definir el número efectivo de componentes k a retener
    if k == "max":
        k = max(1, r)
    else:
        k = min(int(k), len(valores_singulares), r)
        
    # Recortar los valores singulares a k dimensiones
    hatSig = valores_singulares[:k]
    
    # Inicializar las matrices finales truncadas
    hatU = np.zeros((m, k))
    hatV = np.zeros((n, k))
    
    if not transpuesta:
        # Caso m >= n: 
        # res1 es V (n x n), res3 es B = U * Sigma (m x n)
        V = res1
        U_sig = res3
        
        hatV = V[:, :k]
        for i in range(k):
            if hatSig[i].item() > tol:
                hatU[:, i] = U_sig[:, i] / hatSig[i]
            else:
                hatU[:, i] = 0.0
    else:
        # Caso m < n (estaba transpuesta): 
        # res1 es B = V * Sigma (n x m), res3 es U (m x m)
        V_sig = res1
        U = res3
        
        hatU = U[:, :k]
        for i in range(k):
            if hatSig[i] > tol:
                hatV[:, i] = V_sig[:, i] / hatSig[i]
            else:
                hatV[:, i] = 0.0
                
    return hatU, hatSig, hatV

# Matrices al azar
def genera_matriz_para_test(m,n=2,tam_nucleo=0):
    if tam_nucleo == 0:
        A = np.random.random((m,n))
    else:
        A = np.random.random((m,tam_nucleo))
        A = np.hstack([A,A])
    return(A)

def test_svd_reducida_mn(A,tol=1e-15):
    m,n = A.shape
    hU,hS,hV = svd_reducida(A,tol=tol)
    nU,nS,nVT = np.linalg.svd(A)
    r = len(hS)+1
    assert np.all(np.abs(np.abs(np.diag(hU.T @ nU))-1)<10**r*tol), 'Revisar calculo de hat U en ' + str((m,n))
    assert np.all(np.abs(np.abs(np.diag(nVT @ hV))-1)<10**r*tol), 'Revisar calculo de hat V en ' + str((m,n))
    assert len(hS) == len(nS[np.abs(nS)>tol]), 'Hay cantidades distintas de valores singulares en ' + str((m,n))
    assert np.all(np.abs(hS-nS[np.abs(nS)>tol])<10**r*tol), 'Hay diferencias en los valores singulares en ' + str((m,n))

for m in [2,5,10,20]:
    for n in [2,5,10,20]:
        for _ in range(10):
            A = genera_matriz_para_test(m,n)
            test_svd_reducida_mn(A)


# Matrices con nucleo

m = 12
for tam_nucleo in [2,4,6]:
    for _ in range(10):
        A = genera_matriz_para_test(m,tam_nucleo=tam_nucleo)
        test_svd_reducida_mn(A)

# Tamaños de las reducidas
A = np.random.random((8,6))
for k in [1,3,5]:
    hU,hS,hV = svd_reducida(A,k=k)
    assert hU.shape[0] == A.shape[0], 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[0] == A.shape[1], 'Dimensiones de hV incorrectas(caso a)'
    assert hU.shape[1] <= k, 'Dimensiones de hU incorrectas (caso a)'
    assert hV.shape[1] <= k, 'Dimensiones de hV incorrectas(caso a)'
    assert len(hS) <= k, 'Tamaño de hS incorrecto'
