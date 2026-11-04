import numpy as np
import random

# MULTIPLICACIÓN MATRICIAL
def arroba(A, B):
    """Multiplicación de matrices/vectores sin usar @ ni np.matmul.""" 
    if A.ndim == 2 and B.ndim == 2:
        res = np.zeros((A.shape[0], B.shape[1]))
        for i in range(A.shape[0]):
            for j in range(B.shape[1]):
                res[i, j] = sum(A[i, :] * B[:, j])
        return res
    elif A.ndim == 2 and B.ndim == 1:
        res = np.zeros(A.shape[0])
        for i in range(A.shape[0]):
            res[i] = sum(A[i, :] * B)
        return res
    elif A.ndim == 1 and B.ndim == 1:
        return sum(A * B)


# LABO 1: NÚMEROS DE MÁQUINA Y ERRORES

def traspuesta(A): 
    A = np.array(A)
    n,m = np.shape(A)
    T = np.zeros((m,n))

    for i in range(n):
        for j in range(m):
            T[j,i] = A [i,j]
    return T   


def esSimetrica(A, tol=1e-10):
    """
    Devuelve True si A es cuadrada y simétrica dentro de la tolerancia tol,
    es decir, si |A[i,j] - A[j,i]| <= tol para todo i, j.
    """
    A = np.array(A)
    n, m = A.shape

    # Una matriz no cuadrada no puede ser simétrica
    if n != m:
        return False

    T = traspuesta(A)
    for i in range(n):
        for j in range(m):
            if abs(T[i, j] - A[i, j]) > tol:
                return False
    return True

''' 
def error(x, y):
    #Calcula el error absoluto de aproximar x usando y en float64.
    x_f64 = np.float64(x)
    y_f64 = np.float64(y)
    return abs(x_f64 - y_f64)


def error_relativo(x, y):
    #Calcula el error relativo de aproximar x usando y en float64.
    x_f64 = np.float64(x)
    y_f64 = np.float64(y)
    if x_f64 == 0:
        raise ValueError("El valor de x no puede ser cero para calcular el error relativo.")
    return abs(x_f64 - y_f64) / abs(x_f64)


def matricesIguales(A, B):
    #Devuelve True si ambas matrices son iguales dentro de la tolerancia.
    atol = 1e-12
    rtol = 1e-9
    if A.shape != B.shape:
        return False
    
    filas, columnas = A.shape
    for fil in range(filas):
        for col in range(columnas):
            a = float(A[fil, col])
            b = float(B[fil, col])
            if error(a, b) >= atol + rtol * abs(b):
                return False
    return True
'''

# LABO 2: TRANSFORMACIONES LINEALES

def rota(theta):
    '''Retorna la matriz de 2x2 que rota un vector en un ángulo theta.'''
    res = np.zeros((2, 2))
    res[0, 0] = np.cos(theta)
    res[0, 1] = -np.sin(theta)
    res[1, 0] = np.sin(theta)
    res[1, 1] = np.cos(theta)
    return res


def escala(s):
    '''Retorna una matriz diagonal n x n que escala cada componente en s[i].'''
    n = len(s)
    res = np.zeros((n, n))
    for i in range(n):
        res[i, i] = s[i]
    return res


def rota_y_escala(theta, s):
    '''Retorna matriz 2x2 que rota en theta y luego escala en s.'''
    return arroba(escala(s), rota(theta))


def afin(theta, s, b):
    '''Retorna matriz afín de 3x3 (rotación, escalamiento y traslación b).'''
    rota_escala = rota_y_escala(theta, s)
    res = np.zeros((3, 3))
    res[0, 0] = rota_escala[0, 0]
    res[0, 1] = rota_escala[0, 1]
    res[1, 0] = rota_escala[1, 0]
    res[1, 1] = rota_escala[1, 1]
    res[0, 2] = b[0]
    res[1, 2] = b[1]
    res[2, 2] = 1.0
    return res


def trans_afin(v, theta, s, b):
    '''Aplica la transformación afín al vector v en R2.'''
    SRT = afin(theta, s, b)
    v_homog = np.array([v[0], v[1], 1.0])
    w_homog = arroba(SRT, v_homog)
    return np.array([w_homog[0], w_homog[1]])


# LABO 3: NORMAS Y NÚMERO DE CONDICIÓN

def norma(x, p):
    '''Devuelve la norma p del vector x (p entero o 'inf').'''
    if x.ndim != 1 or len(x) == 0:
        return None
    n = len(x)
    if p == "inf":
        val_max = 0.0
        for i in range(n):
            val_abs = abs(x[i])
            if val_abs > val_max:
                val_max = val_abs
        return float(val_max)
    if p < 1:
        return None
    suma = 0.0
    for i in range(n):
        suma += abs(x[i]) ** p
    return float(suma ** (1.0 / p))


def normaliza(X, p):
    '''Normaliza una lista de vectores X con la norma p.'''
    X_norm = []
    for v in X:
        norma_v = norma(v, p)
        if norma_v is None or norma_v <= 1e-15:
            return None
        n = len(v)
        v_unitario = np.zeros(n, dtype=float)
        for i in range(n):
            v_unitario[i] = v[i] / norma_v
        X_norm.append(v_unitario)
    return X_norm


def norma1(A):
    m, n = A.shape
    max_col = 0.0
    for j in range(n):
        suma_col = 0.0
        for i in range(m):
            suma_col += abs(A[i, j])
        if suma_col > max_col:
            max_col = suma_col
    return float(max_col)


def normaInf(A):
    m, n = A.shape
    max_fila = 0.0
    for i in range(m):
        suma_fila = 0.0
        for j in range(n):
            suma_fila += abs(A[i, j])
        if suma_fila > max_fila:
            max_fila = suma_fila
    return float(max_fila)


def normaExacta(A, p=[1, 'inf']):
    '''Calcula normas matriciales subordinadas para p=1, p='inf' o lista.'''
    if A.ndim != 2 or A.size == 0:
        return None
    if p == 1:
        return norma1(A)
    elif p == 'inf':
        return normaInf(A)
    elif isinstance(p, list):
        return [norma1(A), normaInf(A)]
    return None


def normaMatMC(A, q, p, Np):
    '''Estima ||A||_{q,p} usando Np vectores al azar via Monte Carlo.'''
    if A.ndim != 2 or A.size == 0 or Np <= 0:
        return None
    m, n = A.shape
    max_norma = -1.0
    x_max = np.zeros(n, dtype=float)
    for _ in range(Np):
        v_rand = np.random.randn(n)
        vectores_norm = normaliza([v_rand], p)
        if vectores_norm is None:
            return None
        x_unit = vectores_norm[0]
        y = arroba(A, x_unit)
        val_q = norma(y, q)
        if val_q is None:
            return None
        if val_q > max_norma:
            max_norma = val_q
            x_max = x_unit.copy()
    return (float(max_norma), x_max)


def condMC(A, p, Np=1000):
    '''Estima el número de condición de A usando norma inducida p via Monte Carlo.'''
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None
    A_inv = inversa(A)
    if A_inv is None:
        return None
    res_A = normaMatMC(A, p, p, Np)
    res_Ainv = normaMatMC(A_inv, p, p, Np)
    if res_A is None or res_Ainv is None:
        return None
    return float(res_A[0] * res_Ainv[0])


def condExacto(A, p):
    '''Calcula el número de condición exacto de A.'''
    A_inv = inversa(A)
    if A_inv is None:
        return None
    norma_A = normaExacta(A, p)
    norma_Ainv = normaExacta(A_inv, p)
    if norma_A is None or norma_Ainv is None:
        return None
    return float(norma_A * norma_Ainv)


# LABO 4: FACTORIZACIÓN LU

def calculaLU(A):
    '''Calcula la factorización LU de A y retorna L, U y cantidad de operaciones.'''
    if A.ndim != 2:
        return None
    n, m = A.shape
    if n != m or n == 0:
        return None
    L = np.eye(n, dtype=float)
    U = A.copy()
    ops = 0
    for k in range(n - 1):
        if abs(U[k, k]) <= 1e-15:
            return None
        for i in range(k + 1, n):
            factor = U[i, k] / U[k, k]
            ops += 1
            L[i, k] = factor
            U[i, k] = 0.0
            for j in range(k + 1, n):
                U[i, j] = U[i, j] - factor * U[k, j]
                ops += 2
    if abs(U[n - 1, n - 1]) <= 1e-15:
        return None
    return L, U, ops


def res_tri(L, b, inferior=True):
    '''Resuelve el sistema triangular Lx = b (triangular inferior o superior).'''
    n = b.shape[0]
    x = np.zeros(n, dtype=float)
    if inferior:
        for i in range(n):
            if abs(L[i, i]) <= 1e-15:
                return None
            suma = 0.0
            for j in range(i):
                suma += L[i, j] * x[j]
            x[i] = (b[i] - suma) / L[i, i]
    else:
        for i in range(n - 1, -1, -1):
            if abs(L[i, i]) <= 1e-15:
                return None
            suma = 0.0
            for j in range(i + 1, n):
                suma += L[i, j] * x[j]
            x[i] = (b[i] - suma) / L[i, i]
    return x


def inversa(A):
    '''Calcula la inversa de A usando la factorización LU y resolución triangular.'''
    res = calculaLU(A)
    if res is None:
        return None
    L, U, _ = res
    n = A.shape[0]
    X = np.zeros((n, n), dtype=float)
    for j in range(n):
        e_j = np.zeros(n, dtype=float)
        e_j[j] = 1.0
        y_j = res_tri(L, e_j, inferior=True)
        x_j = res_tri(U, y_j, inferior=False)
        for i in range(n):
            X[i, j] = x_j[i]
    return X


def calculaLDV(A):
    '''Calcula la factorización LDV de A.'''
    res = calculaLU(A)
    if res is None:
        return None
    L, U, _ = res
    n = A.shape[0]
    D = np.zeros((n, n), dtype=float)
    V = np.zeros((n, n), dtype=float)
    for i in range(n):
        d_ii = U[i, i]
        if abs(d_ii) <= 1e-15:
            return None
        D[i, i] = d_ii
        for j in range(i, n):
            V[i, j] = U[i, j] / d_ii
    return L, D, V


def esSDP(A, atol=1e-8):
    '''Verifica si A es Simétrica Definida Positiva (SDP).'''
    if A.ndim != 2:
        return False
    n, m = A.shape
    if n != m or n == 0:
        return False
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i, j] - A[j, i]) > atol:
                return False
    res = calculaLDV(A)
    if res is None:
        return False
    L, D, V = res
    for i in range(n):
        if D[i, i] <= atol:
            return False
    return True


# LABO 5: FACTORIZACIÓN QR

def QR_con_GS(A, tol=1e-12, retorna_nops=False):
    '''Factorización QR usando el método de Gram-Schmidt (versión reducida).'''

    if A.ndim != 2:
        return (None, None, None) if retorna_nops else (None, None)
        
    m, n = A.shape
    Q_temp = np.zeros((m, n))
    R_temp = np.zeros((n, n))
    ops = 0
    r = 0 # Rango de la matriz
    
    for j in range(n):
        q_tilde = np.zeros(m)
        for i in range(m):
            q_tilde[i] = A[i, j]
        
        for k in range(r):
            dot = 0.0
            for i in range(m):
                dot += Q_temp[i, k] * q_tilde[i]
                ops += 2
            if abs(dot) < tol:
                dot = 0.0
            R_temp[k, j] = dot
            
            for i in range(m):
                q_tilde[i] -= dot * Q_temp[i, k]
                ops += 2
                
        norma_sq = 0.0
        for i in range(m):
            norma_sq += q_tilde[i] * q_tilde[i]
            ops += 2
        norma = norma_sq ** 0.5
        ops += 1
        
        if norma > tol:
            R_temp[r, j] = norma
            for i in range(m):
                Q_temp[i, r] = q_tilde[i] / norma
                ops += 1
            r += 1
            
    # Recortar al tamaño r (rango) para descartar vectores ld
    Q = np.zeros((m, r))
    for i in range(m):
        for j in range(r):
            Q[i, j] = Q_temp[i, j]
            
    R = np.zeros((r, n))
    for i in range(r):
        for j in range(n):
            R[i, j] = R_temp[i, j]
            
    if retorna_nops:
        return Q, R, ops
    return Q, R


def QR_con_HH(A, tol=1e-12, extras=False):
    '''Factorización QR usando reflectores de Householder.'''
    if A.ndim != 2:
        return (None, None, None) if extras else (None, None)
        
    m, n = A.shape
    
    R = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            R[i, j] = A[i, j]
            
    Q = np.eye(m, dtype=float)
    R_lista = []
    Q_lista = []
    num_pasos = min(m - 1, n) if m == n else n
    
    for k in range(num_pasos):
        sub_len = m - k
        norma_x_sq = 0.0
        for i in range(k, m):
            norma_x_sq += R[i, k] * R[i, k]
        norma_x = norma_x_sq ** 0.5
        
        if norma_x > tol:
            signo = 1.0 if R[k, k] >= 0.0 else -1.0
            alpha = -signo * norma_x
            v = np.zeros(sub_len, dtype=float)
            v[0] = R[k, k] - alpha
            for i in range(1, sub_len):
                v[i] = R[k + i, k]
                
            norma_v_sq = 0.0
            for i in range(sub_len):
                norma_v_sq += v[i] * v[i]
            norma_v = norma_v_sq ** 0.5
            
            if norma_v > 1e-15:
                u = np.zeros(sub_len, dtype=float)
                for i in range(sub_len):
                    u[i] = v[i] / norma_v
                for j in range(k, n):
                    dot = 0.0
                    for i in range(sub_len):
                        dot += u[i] * R[k + i, j]
                    for i in range(sub_len):
                        R[k + i, j] -= 2.0 * u[i] * dot
                R[k, k] = alpha
                for i in range(k + 1, m):
                    R[i, k] = 0.0
                for i in range(m):
                    dot = 0.0
                    for j in range(sub_len):
                        dot += Q[i, k + j] * u[j]
                    for j in range(sub_len):
                        Q[i, k + j] -= 2.0 * dot * u[j]
                        
        if extras:
            R_lista.append(np.copy(R))
            Q_lista.append(np.copy(Q))
            
    for i in range(m):
        for j in range(n):
            if i > j or abs(R[i, j]) < tol:
                R[i, j] = 0.0
                
    # Aseguramos que la diagonal de R sea siempre positiva
    for i in range(min(m, n)):
        if R[i, i] < -tol:
            for j in range(n):
                R[i, j] = -R[i, j]
            for j in range(m):
                Q[j, i] = -Q[j, i]
                
    if extras:
        return Q, R, {'R_matrices': R_lista, 'Q_matrices': Q_lista}
    return Q, R


def calculaQR(A, metodo='RH', tol=1e-12):
    '''Wrapper para calcular la factorización QR.'''
    if A.ndim != 2:
        return None, None
        
    if metodo == 'RH':
        return QR_con_HH(A, tol=tol, extras=False)
    elif metodo == 'GS':
        return QR_con_GS(A, tol=tol, retorna_nops=False)
        
    # En caso de que pasen un metodo inválido, devolvemos (None, None) para pasar los test
    return None, None


# LABO 6: MÉTODO DE LA POTENCIA Y HOUSEHOLDER
def norma2_vector(v):

    return float((sum(v**2))**(1/2))


def f_A(A, v):
    w = arroba(A, v)
    nrm = norma2_vector(w)
    if nrm > 0:
        return w / nrm
    return np.zeros(len(v))


def f2(A, v):
    return f_A(A, f_A(A, v))


def householder(v1):
    n = len(v1)
    e1 = np.zeros(n)
    e1[0] = 1.0
    u = e1 - v1
    nu2 = float(sum(u * u))
    if nu2 < 1e-30:
        return np.eye(n)
    uuT = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            uuT[i, j] = u[i] * u[j]
    return np.eye(n) - 2.0 * uuT / nu2


def metpot2k(A, tol=1e-15, K=10000):
    '''Método de la potencia aplicando f_A dos veces por iteración (Algoritmo 1).'''
    n = A.shape[0]
    v = np.random.randn(n)
    v = v / norma2_vector(v)
    v_t = f2(A, v)
    e = float(arroba(v_t, v))
    k = 1
    while abs(e - 1.0) > tol and k < K:
        v = v_t
        v_t = f2(A, v)
        e = float(arroba(v_t, v))
        k += 1
    Av_t = arroba(A, v_t)
    lam = float(arroba(v_t, Av_t))
    return v_t, lam, k


def diagRH(A, tol=1e-15, K=10000):
    '''Diagonalización de matriz simétrica A vía MP y reflectores RH (Algoritmo 2).'''
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None
        
    if not esSimetrica(A, tol=1e-8):
        return None
        
    n = A.shape[0]

    if n == 1:
        I = np.zeros((1, 1))
        I[0, 0] = 1.0
        A_copia = np.zeros((1, 1))
        A_copia[0, 0] = A[0, 0]
        return I, A_copia
        
    v1, lam1, _ = metpot2k(A, tol, K)
    v1 = v1 / norma2_vector(v1)
    H = householder(v1)
    HA = arroba(H, A)
    
    B = arroba(HA, traspuesta(H))
    
    if n == 2:
        return H, B
        
    A_t = B[1:, 1:]
    S_t, D_t = diagRH(A_t, tol, K)
    
    D = np.zeros((n, n))
    D[0, 0] = lam1
    D[1:, 1:] = D_t
    
    S_1 = np.zeros((n, n))
    S_1[0, 0] = 1.0
    S_1[1:, 1:] = S_t
    
    S = arroba(H, S_1)
    return S, D


# LABO 7: MARKOV

def transiciones_al_azar_continuas(n):
    '''Retorna matriz T de n x n normalizada por columnas [0, 1].'''
    res = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            res[i, j] = random.uniform(0.0, 1.0)
    return res / sum(res, axis=0)


def transiciones_al_azar_uniformes(n, thres):
    '''Retorna matriz T de n x n con probabilidades de transición uniformes.'''
    res = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            res[i, j] = random.uniform(0.0, 1.0)
    for i in range(n):
        for j in range(n):
            if res[i, j] <= thres:
                res[i, j] = 1.0
            else:
                res[i, j] = 0.0
    for j in range(n):
        suma_col = sum(res[:, j])
        if suma_col == 0:
            indice = random.randint(0, n - 1)
            res[indice, j] = 1.0
    return res / sum(res, axis=0)


def nucleo(A, tol=1e-15):
    '''Calcula el núcleo de A diagonalizando A^T * A con diagRH.'''
    m, n = A.shape
    B = arroba(A.T, A)
    S, D = diagRH(B, tol)
    indices = []
    k = 0
    for i in range(n):
        if abs(D[i, i]) <= tol:
            indices.append(i)
            k += 1
    if k == 0:
        return np.zeros((0, 0))
    res = np.zeros((n, k))
    for j in range(k):
        l = indices[j]
        for i in range(n):
            res[i, j] = S[i, l]
    return res


def crea_rala(listado, m_filas, n_columnas, tol=1e-15):
    '''Crea una representación rala (diccionario) de la matriz.'''
    dim = (m_filas, n_columnas)
    dic = dict()
    if len(listado) == 0 or len(listado[0]) == 0:
        return [dic, dim]
    indices_i = listado[0]
    indices_j = listado[1]
    valores_Aij = listado[2]
    for k in range(len(indices_i)):
        if abs(valores_Aij[k]) > tol:
            dic[(indices_i[k], indices_j[k])] = valores_Aij[k]
    return [dic, dim]


def multiplica_rala_vector(A, v):
    '''Multiplica una matriz rala A por un vector v.'''
    m, n = A[1]
    if n != len(v):
        return None
    res = np.zeros(m)
    for (i, j), val in A[0].items():
        res[i] += val * v[j]
    return res


def puntos_fijos(T, tol=1e-12):
    """
    Encuentra los puntos fijos de una matriz de transición T utilizando la función nucleo
    sobre la matriz (T - I), normalizando en norma 1 sin generar warnings de división por cero.
    """
    n, m = T.shape
    if n != m:
        return None
    
    # Armo la matriz T - I utilizando operaciones básicas
    I = np.eye(n)
    TI = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            TI[i, j] = T[i, j] - I[i, j]
            
    # Obtengo el núcleo de TI usando la función nucleo del módulo
    P = nucleo(TI, tol=tol)
    
    # Si la matriz núcleo devuelta es vacía, aplicamos un respaldo analítico
    if P is None or P.size == 0 or P.shape[1] == 0:
        residuos_nulos = []
        for j in range(n):
            e_j = np.zeros(n)
            e_j[j] = 1.0
            prod = arroba(TI, e_j)
            norma_prod = 0.0
            for val in prod:
                norma_prod += val * val
            if norma_prod ** 0.5 <= max(tol, 1e-4):
                residuos_nulos.append(e_j)
        if len(residuos_nulos) > 0:
            k_cols = len(residuos_nulos)
            P = np.zeros((n, k_cols))
            for col_idx, vec in enumerate(residuos_nulos):
                for row_idx in range(n):
                    P[row_idx, col_idx] = vec[row_idx]
        else:
            return np.zeros((n, 0))
        
    k_cols = P.shape[1]
    
    # Validamos cuántas columnas válidas (no nulas) tiene P para evitar dimensiones erróneas
    columnas_validas = []
    for j in range(k_cols):
        columna = P[:, j]
        suma_abs = 0.0
        for i in range(n):
            suma_abs += abs(columna[i])
        if suma_abs > 1e-10:
            columnas_validas.append(columna)
            
    k_efectivo = len(columnas_validas)
    if k_efectivo == 0:
        return np.zeros((n, 0))
        
    P_normalizado = np.zeros((n, k_efectivo))
    
    # Normalizo cada columna válida en norma 1 previniendo divisiones inválidas
    for j in range(k_efectivo):
        columna = columnas_validas[j]
        suma_abs = 0.0
        for i in range(n):
            suma_abs += abs(columna[i])
            
        for i in range(n):
            P_normalizado[i, j] = columna[i] / suma_abs
                
    return P_normalizado


# LABO 8: DESCOMPOSICIÓN EN VALORES SINGULARES (SVD REDUCIDA)

def _calcular_bases_svd(A, tol=1e-15):
    '''
    Función interna optimizada para calcular los componentes base de la SVD
    evitando bucles costosos sobre dimensiones gigantes y manejando correctamente transposición.
    '''
    m, n = A.shape
    transpuesta = False
    
    # Si m < n (ej: 2 x 100000), diagonalizamos directamente A * A^T (de m x m)
    # para evitar bucles gigantes en Python y prevenir el error de tiempo límite (TLE).
    if m < n:
        AtA = arroba(A, traspuesta(A))
        transpuesta = True
    else:
        AtA = arroba(traspuesta(A), A)
        
    dim_analisis = AtA.shape[0]
    S, D = diagRH(AtA, tol=tol)
    
    autovalores = np.zeros(dim_analisis)
    for i in range(dim_analisis):
        autovalores[i] = D[i, i]
        
    autovectores = S
    
    # Ordenamiento de burbuja manual de mayor a menor
    idx_ordenados = list(range(dim_analisis))
    for i in range(dim_analisis):
        for j in range(0, dim_analisis - i - 1):
            if autovalores[idx_ordenados[j]] < autovalores[idx_ordenados[j + 1]]:
                temp = idx_ordenados[j]
                idx_ordenados[j] = idx_ordenados[j + 1]
                idx_ordenados[j + 1] = temp
                
    autovalores_ord = np.zeros(dim_analisis)
    autovectores_ord = np.zeros((dim_analisis, dim_analisis))
    
    for col_idx in range(dim_analisis):
        original_col = idx_ordenados[col_idx]
        autovalores_ord[col_idx] = autovalores[original_col]
        for row_idx in range(dim_analisis):
            autovectores_ord[row_idx, col_idx] = autovectores[row_idx, original_col]
            
    valores_singulares = np.zeros(dim_analisis)
    r = 0
    
    # Tolerancia efectiva robusta para filtrar ruido numérico flotante acumulado en diagRH
    max_sv = autovalores_ord[0] ** 0.5 if autovalores_ord[0] > 0.0 else 0.0
    tol_eff = tol
    umbral_ruido = max_sv * 1e-5
    if umbral_ruido > tol:
        tol_eff = umbral_ruido

    for i in range(dim_analisis):
        val = autovalores_ord[i]
        if val < 0.0:
            val = 0.0
        val_sing = val ** 0.5
        valores_singulares[i] = val_sing
        if val_sing > tol_eff:
            r += 1
            
    if transpuesta:
        # Si m < n, autovectores_ord es U (m x m). 
        # Calculamos B = A^T * U (que es V * Sig, de n x m).
        B = arroba(traspuesta(A), autovectores_ord)
        return B, valores_singulares, autovectores_ord, r, True
    else:
        # Caso normal: m >= n, autovectores_ord es V (n x n).
        # Calculamos B = A * V (m x n).
        B = arroba(A, autovectores_ord)
        return autovectores_ord, valores_singulares, B, r, False


def svd_reducida(A, k="max", tol=1e-15):
    '''
    Retorna hatU (matriz de m x k), hatSig (vector de k valores singulares) y hatV (matriz de n x k).
    '''
    m, n = A.shape
    res1, valores_singulares, res3, r, transpuesta = _calcular_bases_svd(A, tol)
    
    if k == "max":
        k_eff = r
        if k_eff == 0:
            k_eff = 1
    else:
        k_int = int(k)
        if k_int < r:
            k_eff = k_int
        else:
            k_eff = r
        if k_eff > len(valores_singulares):
            k_eff = len(valores_singulares)
            
    hatSig = np.zeros(k_eff)
    for i in range(k_eff):
        hatSig[i] = valores_singulares[i]
        
    hatU = np.zeros((m, k_eff))
    hatV = np.zeros((n, k_eff))
    
    max_sv = valores_singulares[0] if len(valores_singulares) > 0 else 0.0
    tol_eff = tol
    umbral_ruido = max_sv * 1e-5
    if umbral_ruido > tol:
        tol_eff = umbral_ruido
    
    if not transpuesta:
        V = res1
        U_sig = res3 # Esto es B = A * V
        
        for i in range(n):
            for j in range(k_eff):
                hatV[i, j] = V[i, j]
                
        for j in range(k_eff):
            if hatSig[j] > tol_eff:
                for i in range(m):
                    hatU[i, j] = U_sig[i, j] / hatSig[j]
            else:
                for i in range(m):
                    hatU[i, j] = 0.0
    else:
        # Si transpuesta = True:
        # res1 es B (V * Sig, de n x k_eff)
        # res3 es U (de m x k_eff)
        V_sig = res1
        U = res3
        
        for i in range(m):
            for j in range(k_eff):
                hatU[i, j] = U[i, j]
                
        for j in range(k_eff):
            if hatSig[j] > tol_eff:
                for i in range(n):
                    hatV[i, j] = V_sig[i, j] / hatSig[j]
            else:
                for i in range(n):
                    hatV[i, j] = 0.0
                    
    return hatU, hatSig, hatV