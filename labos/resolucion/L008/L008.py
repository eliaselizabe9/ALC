import numpy as np



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


def svd_reducida(id_matriz, k="max", tol=1e-15):
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
            if hatSig[i] > tol:
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