import numpy as np

# ==============================================================================
# FUNCIONES AUXILIARES L006 (metpot2k y diagRH)
# ==============================================================================

"""
# Borrador inicial metpot2k: intento iterativo acumulando potencias directas
def metpot2k(A, tol=1e-15, K=1000):
    n = A.shape[0]
    v = np.random.rand(n)
    v = v / np.sqrt(v @ v)
    
    # Intento aplicando A^2 directamente en un solo paso
    A2 = np.matmul(A, A)
    v_tilde = A2 @ v
    v_tilde = v_tilde / np.sqrt(v_tilde @ v_tilde)
    
    k = 0
    while k < K:
        v_next = A2 @ v_tilde
        norma = np.sqrt(v_next @ v_next)
        if norma < 1e-15:
            break
        v_next = v_next / norma
        if abs(abs(v_next @ v_tilde) - 1.0) < tol:
            v_tilde = v_next
            break
        v_tilde = v_next
        k += 1
        
    l = v_tilde @ (A @ v_tilde)
    return v_tilde, l, k
"""

def metpot2k(A, tol=1e-15, K=1000):
    n = A.shape[0]
    v = np.random.rand(n)
    norm_v = np.sqrt(v @ v)
    if norm_v > 0:
        v = v / norm_v

    v_tilde = A @ (A @ v)
    norm_vt = np.sqrt(v_tilde @ v_tilde)
    if norm_vt > 0:
        v_tilde = v_tilde / norm_vt
    else:
        return v, 0.0, 0

    e = v_tilde @ v
    k = 0
    while abs(e - 1.0) > tol and k < K:
        v = v_tilde
        v_tilde = A @ (A @ v)
        norm_vt = np.sqrt(v_tilde @ v_tilde)
        if norm_vt > 0:
            v_tilde = v_tilde / norm_vt
        else:
            return v, 0.0, k
        e = v_tilde @ v
        k += 1

    Av = A @ v_tilde
    l = v_tilde @ Av
    return v_tilde, l, k


"""
# Borrador Householder intentando proyectar con signo directo de v1[0]
def diagRH(A, tol=1e-15, K=1000):
    n = A.shape[0]
    if n == 1:
        return np.eye(1), A.copy()
    
    v1, l1, _ = metpot2k(A, tol, K)
    e1 = np.zeros(n)
    e1[0] = 1.0
    
    # Alternativa sin signo adaptativo
    u = v1 - e1
    nu = np.sqrt(u @ u)
    if nu > 1e-15:
        u = u / nu
        H = np.eye(n) - 2.0 * np.outer(u, u)
    else:
        H = np.eye(n)
        
    B = H @ A @ H
    S_rec, D_rec = diagRH(B[1:, 1:], tol, K)
    
    D = np.zeros((n, n))
    D[0, 0] = l1
    D[1:, 1:] = D_rec
    
    # Acumulacion directa de bases
    S = np.eye(n)
    S[1:, 1:] = S_rec
    S = np.matmul(H, S)
    return S, D
"""

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


# ==============================================================================
# FUNCIONES PRINCIPALES
# ==============================================================================

"""
# Borrador transiciones_azar_uniforme normalizando por fila
def transiciones_azar_uniforme(n, thres=1.0):
    R = np.random.uniform(0.0, 1.0, (n, n))
    Mascara = (R <= thres).astype(float)
    T = R * Mascara
    sumas = np.matmul(T, np.ones((n, 1)))
    for i in range(n):
        if sumas[i, 0] > 1e-15:
            T[i, :] = T[i, :] / sumas[i, 0]
        else:
            T[i, :] = 1.0 / n
    return T
"""

def transiciones_azar_uniforme(n, thres=1.0):
    R = np.random.uniform(0.0, 1.0, (n, n))
    Mascara = (R <= thres).astype(float)
    T = R * Mascara
    
    # Suma por columnas calculada mediante producto vectorial: e.T @ T
    sumas_col = np.ones(n) @ T
    
    for j in range(n):
        if sumas_col[j] > 1e-15:
            T[:, j] = T[:, j] / sumas_col[j]
        else:
            T[:, j] = 1.0 / n
    return T


def transiciones_al_azar_continuas(n):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    Retorna matriz T de n x n normalizada por columnas, y con entradas al azar en el intervalo [0,1]
    """
    return transiciones_azar_uniforme(n, thres=1.0)


"""
# Borrador transiciones_azar acumulando mascara booleana directa
def transiciones_azar(n, thres):
    R = np.random.uniform(0.0, 1.0, (n, n))
    T = (R <= thres).astype(float)
    k_counts = np.ones(n) @ T
    for j in range(n):
        if k_counts[j] > 0:
            T[:, j] = T[:, j] * (1.0 / k_counts[j])
        else:
            T[:, j] = 0.0
    return T
"""

def transiciones_azar(n, thres):
    R = np.random.uniform(0.0, 1.0, (n, n))
    T = (R <= thres).astype(float)
    
    # Conteo de no nulos por columna via e.T @ T
    k_counts = np.ones(n) @ T
    for j in range(n):
        if k_counts[j] > 0:
            T[:, j] = T[:, j] / k_counts[j]
        else:
            T[:, j] = 1.0 / n
    return T


def transiciones_al_azar_uniformes(n, thres):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    thres probabilidad de que una entrada sea distinta de cero.
    Retorna matriz T de n x n normalizada por columnas.
    """
    return transiciones_azar(n, thres)


"""
# Borrador nucleo resolviendo autovalores sobre A @ A.T
def nucleo(A, tol=1e-15):
    m, n = A.shape
    M = np.matmul(A, A.T)
    res = diagRH(M, tol=tol, K=1000)
    if res is None:
        raise ValueError("Error en la diagonalizacion")
    S, D = res
    idx = [i for i in range(m) if abs(D[i, i]) <= tol]
    if len(idx) == 0:
        return np.zeros((0, 0))
    return S[:, idx]
"""

def nucleo(A, tol=1e-15):
    """
    A una matriz de m x n
    Calcula el nucleo de la matriz A diagonalizando A.T @ A con diagRH.
    """
    # M = A^T * A mediante np.matmul
    M = np.matmul(A.T, A)

    res = diagRH(M, tol=tol, K=1000)
    if res is None:
        raise ValueError("Falla en la diagonalizacion de A^T * A")

    S, D = res

    indices_nucleo = [j for j in range(A.shape[1]) if abs(D[j, j]) <= tol]

    k = len(indices_nucleo)
    if k == 0:
        return np.zeros((0, 0))

    return S[:, indices_nucleo]


"""
# Borrador crea_rala parseando como tuplas desordenadas sin filtro de tolerancia
def crea_rala(listado, m_filas, n_columnas, tol=1e-15):
    dic_rala = {}
    if len(listado) == 3 and len(listado[0]) > 0:
        for idx in range(len(listado[2])):
            i = int(listado[0][idx])
            j = int(listado[1][idx])
            dic_rala[(i, j)] = listado[2][idx]
    return [dic_rala, (m_filas, n_columnas)]
"""

def crea_rala(listado, m_filas, n_columnas, tol=1e-15):
    dic_rala = {}
    if len(listado) == 3:
        I = listado[0]
        J = listado[1]
        V = listado[2]
        num_elementos = len(V)

        for k in range(num_elementos):
            val = V[k]
            if abs(val) >= tol:
                i = int(I[k])
                j = int(J[k])
                dic_rala[(i, j)] = val

    dimensiones = (m_filas, n_columnas)
    return [dic_rala, dimensiones]


"""
# Borrador multiplica_rala_vector iterando con dot_vec producto interno por fila
def multiplica_rala_vector(A, v):
    dic_rala, (m_filas, n_columnas) = A
    if len(v) != n_columnas:
        raise ValueError("Dimensiones incompatibles")
    w = np.zeros(m_filas)
    for fil in range(m_filas):
        fila_vec = np.zeros(n_columnas)
        for col in range(n_columnas):
            if (fil, col) in dic_rala:
                fila_vec[col] = dic_rala[(fil, col)]
        w[fil] = fila_vec @ v
    return w
"""

def multiplica_rala_vector(A, v):
    dic_rala, (m_filas, n_columnas) = A

    if len(v) != n_columnas:
        raise ValueError("Dimensiones incompatibles entre matriz rala y vector")

    # Reconstrucción densa temporal para operar matricialmente con @ y np.matmul
    A_densa = np.zeros((m_filas, n_columnas))
    for (i, j), a_ij in dic_rala.items():
        A_densa[i, j] = a_ij

    return np.matmul(A_densa, v)


# ==============================================================================
# FUNCIONES DE VERIFICACIÓN Y TESTS
# ==============================================================================

def es_markov(T, tol=1e-6):
    n = T.shape[0]
    for i in range(n):
        for j in range(n):
            if T[i, j] < 0:
                return False
    # Suma de columnas verificada mediante el producto 1.T @ T
    sumas_col = np.ones(n) @ T
    for j in range(n):
        if np.abs(sumas_col[j] - 1.0) > tol:
            return False
    return True


def es_markov_uniforme(T, thres=1e-6):
    if not es_markov(T, thres):
        return False
    m = T.shape[1]
    for j in range(m):
        non_zero = T[:, j][T[:, j] > thres]
        close = all(np.abs(non_zero - non_zero[0]) < thres)
        if not close:
            return False
    return True


def esNucleo(A, S, tol=1e-5):
    for col in S.T:
        res = A @ col
        if not np.allclose(res, np.zeros(A.shape[0]), atol=tol):
            return False
    return True


if __name__ == "__main__":
    for i in range(1, 100):
        T = transiciones_al_azar_continuas(i)
        assert es_markov(T), f"transiciones_al_azar_continuas fallo para n={i}"
        
        T = transiciones_al_azar_uniformes(i, 0.3)
        assert es_markov_uniforme(T), f"transiciones_al_azar_uniformes fallo para n={i}"
        T = transiciones_al_azar_uniformes(i, 0.01)
        assert es_markov_uniforme(T), f"transiciones_al_azar_uniformes fallo para n={i}"
        
    A = np.eye(3)
    S = nucleo(A)
    assert S.shape[0] == 0, "nucleo fallo para matriz identidad"
    A[1, 1] = 0
    S = nucleo(A)
    msg = "nucleo fallo para matriz con un cero en diagonal"
    assert esNucleo(A, S), msg
    assert S.shape == (3, 1), msg
    assert abs(S[2, 0]) < 1e-2, msg
    assert abs(S[0, 0]) < 1e-2, msg

    v = np.random.random(5)
    v = v / np.sqrt(v @ v)
    H = np.eye(5) - np.outer(v, v)
    S = nucleo(H)
    msg = "nucleo fallo para matriz de proyeccion ortogonal"
    assert S.shape == (5, 1), msg
    v_gen = S[:, 0]
    v_gen = v_gen / np.sqrt(v_gen @ v_gen)
    assert np.allclose(v, v_gen) or np.allclose(v, -v_gen), msg

    listado = [[0, 17], [3, 4], [0.5, 0.25]]
    A_rala_dict, dims = crea_rala(listado, 32, 89)
    assert dims == (32, 89)
    assert A_rala_dict[(0, 3)] == 0.5
    assert A_rala_dict[(17, 4)] == 0.25
    assert len(A_rala_dict) == 2

    listado = [[32, 16, 5], [3, 4, 7], [7, 0.5, 0.25]]
    A_rala_dict, dims = crea_rala(listado, 50, 50)
    assert dims == (50, 50)
    assert A_rala_dict.get((32, 3)) == 7
    assert A_rala_dict[(16, 4)] == 0.5
    assert A_rala_dict[(5, 7)] == 0.25

    listado = [[1, 2, 3], [4, 5, 6], [1e-20, 0.5, 0.25]]
    A_rala_dict, dims = crea_rala(listado, 10, 10)
    assert dims == (10, 10)
    assert (1, 4) not in A_rala_dict
    assert A_rala_dict[(2, 5)] == 0.5
    assert A_rala_dict[(3, 6)] == 0.25
    assert len(A_rala_dict) == 2

    listado = []
    A_rala_dict, dims = crea_rala(listado, 10, 10)
    assert dims == (10, 10)
    assert len(A_rala_dict) == 0

    listado = [[0, 1, 2], [0, 1, 2], [1, 2, 3]]
    A_rala = crea_rala(listado, 3, 3)
    v = np.random.random(3)
    v = v / np.sqrt(v @ v)
    res = multiplica_rala_vector(A_rala, v)
    A = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    res_esperado = A @ v
    assert np.allclose(res, res_esperado)

    A = np.random.random((5, 5))
    A = A * (A > 0.5) 
    listado = [[], [], []]
    for i in range(5):
        for j in range(5):
            listado[0].append(i)
            listado[1].append(j)
            listado[2].append(A[i, j])
            
    A_rala = crea_rala(listado, 5, 5)
    v = np.random.random(5)
    assert np.allclose(multiplica_rala_vector(A_rala, v), A @ v)
    print("¡Todos los tests pasaron exitosamente!")