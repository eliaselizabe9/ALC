import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#import random

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


def error(x, y):
    #Calcula el error absoluto de aproximar x usando y en float64.
    x_f64 = float(x)
    y_f64 = float(y)
    return abs(x_f64 - y_f64)


def error_relativo(x, y):
    #Calcula el error relativo de aproximar x usando y en float64.
    x_f64 = float(x)
    y_f64 = float(y)
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
    U = np.array(A, dtype=float)
    
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
    # Protegemos contra divisiones por números extremadamente pequeños
    if nrm > 1e-15:
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
    
    # Inicialización determinística asimétrica para evitar caer en el núcleo 
    # de matrices de suma cero (como las de Markov).
    v = np.zeros(n, dtype=float)
    for i in range(n):
        v[i] = float(i + 1)
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
    
    # para evitar la división por ceros:
    nrm_v1 = norma2_vector(v1)
    if nrm_v1 > 1e-15:
        v1 = v1 / nrm_v1
    else:
        # Si la matriz era de ceros, el vector es nulo. 
        # Cualquier vector es autovector, usamos el canónico e1.
        v1 = np.zeros(n)
        v1[0] = 1.0
    
    H = householder(v1)
    HA = arroba(H, A)
    
    B = arroba(HA, traspuesta(H))
    B = 0.5 * (B + traspuesta(B)) # Cancelado de error residual de redondeo de float ~ 10^-7
    
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
    '''Retorna matriz T de n x n normalizada por columnas en [0, 1].'''
    res = np.zeros((n, n))
    
    # 1. Generar elementos al azar en [0, 1] usando NumPy
    for i in range(n):
        for j in range(n):
            res[i, j] = np.random.uniform(0.0, 1.0)
            
    # 2. Normalización manual por columnas sin utilizar axis
    for j in range(n):
        suma_col = 0.0
        for i in range(n):
            suma_col += res[i, j]
            
        if suma_col > 1e-15:
            for i in range(n):
                res[i, j] = res[i, j] / suma_col
                
    return res

def transiciones_al_azar_uniformes(n, thres):
    '''Retorna matriz T de n x n con probabilidades de transición uniformes.'''
    res = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            val = np.random.uniform(0.0, 1.0)
            if val <= thres:
                res[i, j] = 1.0
            else:
                res[i, j] = 0.0
                
    # Normalización manual por columnas
    for j in range(n):
        suma_col = 0.0
        for i in range(n):
            suma_col += res[i, j]
            
        # Si una columna quedó completamente en ceros, asignamos un 1 a una fila al azar
        if suma_col == 0.0:
            idx = np.random.randint(0, n)
            res[idx, j] = 1.0
            suma_col = 1.0
            
        for i in range(n):
            res[i, j] = res[i, j] / suma_col
            
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

'''
def nucleo(A, tol=1e-15):
    #Calcula el núcleo de A diagonalizando A^T * A con diagRH.
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
'''
def nucleo(A, tol=1e-15):
    '''Calcula el núcleo de A diagonalizando A^T * A con diagRH.'''
    m, n = A.shape
    B = arroba(traspuesta(A), A)
    S, D = diagRH(B, tol)
    indices = []
    k = 0
    for i in range(n):
        val_aut = D[i, i]
        if val_aut < 0.0:
            val_aut = 0.0
        if val_aut <= tol:
            indices.append(i)
            k += 1
    if k == 0:
        return np.zeros((n, 0))
    res = np.zeros((n, k))
    for j in range(k):
        l = indices[j]
        for i in range(n):
            res[i, j] = S[i, l]
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
    Función interna optimizada para calcular los componentes base de la SVD.
    Se utiliza vectorización 1D por columnas para la multiplicación final,
    evitando el error de Tiempo Excedido (TLE) en matrices de 100000 componentes.
    '''
    m, n = A.shape
    transpuesta = False
    
    # Si m < n, diagonalizamos A * A^T
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
        # Asignación de columnas vectorizada para velocidad máxima
        autovectores_ord[:, col_idx] = autovectores[:, original_col]
            
    valores_singulares = np.zeros(dim_analisis)
    r = 0
    
    # Tolerancia efectiva robusta para filtrar ruido numérico
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
        # B = A^T * autovectores_ord
        # Vectorizamos la multiplicación en 1D para que procesar 100000 iteraciones tarde milisegundos
        T_A = traspuesta(A)
        k = autovectores_ord.shape[1]
        B = np.zeros((T_A.shape[0], k))
        for j in range(k):
            col_res = np.zeros(T_A.shape[0])
            for p in range(T_A.shape[1]):
                col_res += T_A[:, p] * autovectores_ord[p, j]
            B[:, j] = col_res
            
        return B, valores_singulares, autovectores_ord, r, True
    else:
        # B = A * autovectores_ord
        k = autovectores_ord.shape[1]
        B = np.zeros((A.shape[0], k))
        for j in range(k):
            col_res = np.zeros(A.shape[0])
            for p in range(A.shape[1]):
                col_res += A[:, p] * autovectores_ord[p, j]
            B[:, j] = col_res
            
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
        U_sig = res3 # B = A * V
        
        # Asignación en bloque rápido
        hatV[:, :k_eff] = V[:, :k_eff]
        
        # División por vectores completos para máxima eficiencia
        for j in range(k_eff):
            if hatSig[j] > tol_eff:
                hatU[:, j] = U_sig[:, j] / hatSig[j]
            else:
                hatU[:, j] = 0.0
    else:
        # transpuesta = True:
        V_sig = res1 # B = A^T * U
        U = res3
        
        hatU[:, :k_eff] = U[:, :k_eff]
        
        for j in range(k_eff):
            if hatSig[j] > tol_eff:
                hatV[:, j] = V_sig[:, j] / hatSig[j]
            else:
                hatV[:, j] = 0.0
                    
    return hatU, hatSig, hatV

# FUNCIONES EJ 2

# EJERCICIO 2a
def obtenerQ3d(imgPts, objPts):
    """
    Construye la matriz Q de dimensión 2n x 12 para el sistema homogéneo Qm = 0.
    
    Parámetros:
        imgPts: array de n x 2 con coordenadas proyectadas (u, v) en la imagen.
        objPts: array de n x 3 con coordenadas reales (U, V, W) en el espacio 3D.
    Retorna:
        Q: matriz de 2n x 12 con las restricciones apiladas.
    """
    n = objPts.shape[0]
    Q = np.zeros((2 * n, 12))
    
    for i in range(n):
        U, V, W = objPts[i, 0], objPts[i, 1], objPts[i, 2]
        u, v = imgPts[i, 0], imgPts[i, 1]
        
        # Armo las coordenadas homogéneas del punto 3D
        P_i = np.array([U, V, W, 1.0])
        
        # Armo la fila 2i correspondiente a la restricción en la coordenada u
        Q[2 * i, 0:4] = P_i
        Q[2 * i, 4:8] = 0.0
        Q[2 * i, 8:12] = -u * P_i
        
        # Armo la fila 2i + 1 correspondiente a la restricción en la coordenada v
        Q[2 * i + 1, 0:4] = 0.0
        Q[2 * i + 1, 4:8] = P_i
        Q[2 * i + 1, 8:12] = -v * P_i
        
    return Q

# EJERCICIO 2b
def resolverDLT3D(Q):
    """
    Resuelvo el sistema DLT 3D buscando el vector m que minimiza ||Qm||^2 sujeto a ||m|| = 1.
    Diagonalizo Q^T * Q usando diagRH de alc y extraigo el autovector del menor autovalor.
    """
    # Calculo la matriz simétrica Q^T * Q de 12x12
    QtQ = arroba(Q.T, Q)
    
    # Diagonalizo QtQ mediante el Método de la Potencia y reflectores de Householder
    S, D = diagRH(QtQ, tol=1e-15, K=10000)
    
    # Busco el índice del autovalor de menor valor absoluto
    idx_min = 0
    min_val = abs(D[0, 0])
    for i in range(1, 12):
        val_abs = abs(D[i, i])
        if val_abs < min_val:
            min_val = val_abs
            idx_min = i
            
    # Extraigo el autovector m asociado y lo normalizo en norma 2
    m = S[:, idx_min]
    norma_m = np.sqrt(np.sum(m**2))
    m = m / norma_m
    
    # Reorganizo el vector m de 12 elementos en la matriz de proyección M de 3x4
    M = m.reshape((3, 4))

    # Fijamos el signo de M para que det(M[:, :3]) > 0
    A0 = M[:, :3]
    det_M3 = (A0[0, 0] * (A0[1, 1] * A0[2, 2] - A0[1, 2] * A0[2, 1]) -
              A0[0, 1] * (A0[1, 0] * A0[2, 2] - A0[1, 2] * A0[2, 0]) +
              A0[0, 2] * (A0[1, 0] * A0[2, 1] - A0[1, 1] * A0[2, 0]))
    if det_M3 < 0:
        M = -M

    return M


from alc import QR_con_HH
# EJERCICIOS 2.c y 2.d: DESCOMPOSICIÓN DE M EN K, R Y t
def descomposicion_RQ_3x3(A):
    """
    Realizo la descomposición RQ de una matriz A de 3x3 mediante la factorización QR de alc
    Descompongo A en K * R, donde K es triangular superior y R es ortogonal.
    """
    # Defino la matriz de permutación/reversión P de 3x3
    P = np.array([[0., 0., 1.],
                  [0., 1., 0.],
                  [1., 0., 0.]])
    
    # Mapeo el problema a una factorización QR mediante (P * A)^T
    PA_T = arroba(P, A).T
    
    # Factorizo utilizando QR_con_HH de alc
    Q_qr, R_qr = QR_con_HH(PA_T)
    
    # Reconstruyo la matriz K (triangular superior) y R (ortogonal)
    K = arroba(P, arroba(R_qr.T, P))
    R = arroba(P, Q_qr.T)
    
    # Ajusto los signos de las columnas de K y filas de R para asegurar diagonal positiva en K
    for i in range(3):
        if K[i, i] < 0:
            K[:, i] = -K[:, i]
            R[i, :] = -R[i, :]
            
    # Calculo el determinante de R para garantizar det(R) = +1 (rotación pura)
    det_R = (R[0, 0] * (R[1, 1] * R[2, 2] - R[1, 2] * R[2, 1]) -
             R[0, 1] * (R[1, 0] * R[2, 2] - R[1, 2] * R[2, 0]) +
             R[0, 2] * (R[1, 0] * R[2, 1] - R[1, 1] * R[2, 0]))
             
    if det_R < 0:
        R = -R
        K = -K
        
    return K, R


def devolverMatricesCalibracion3d(M, fijar_sesgo_cero=False, igualar_focales=False):
    """
    Descompongo la matriz de proyección M (3x4) en las matrices K (intrínsecos), R (rotación)
    y el vector t (traslación).
    Permito opcionalmente imponer la condición de sesgo nulo (gamma = 0)
    y/o igualar las distancias focales (alpha_u = alpha_v).
    """
    # Separo la submatriz A (3x3) y la columna b (3x1) de M = [A | b]
    A = M[:, 0:3]
    b = M[:, 3]
    
    # Descompongo la submatriz A en K y R usando la rutina RQ
    K, R = descomposicion_RQ_3x3(A)

    # Si descomposicion_RQ_3x3 invirtió K y R por det_R < 0, el producto K @ R
    # queda con signo opuesto a A (producto interno negativo) ---> invertimos b también.
    KR = arroba(K, R)
    if np.sum(A * KR) < 0:
        b = -b


    # Hay que dividir por el factor de escala s así que lo defino (es K[2, 2] antes de normalizar)
    s = K[2, 2]
    
    # Normalizo K para asegurar que el elemento homogéneo K[2,2] valga 1
    K = K / s

    # Impongo las restricciones opcionales sobre la matriz de intrínsecos K
    if fijar_sesgo_cero:
        K[0, 1] = 0.0  # gamma = 0
        
    if igualar_focales:
        alpha_prom = (K[0, 0] + K[1, 1]) / 2.0
        K[0, 0] = alpha_prom
        K[1, 1] = alpha_prom

    # Escalamos el vector b con el mismo factor s para recuperar la unidad métrica real de t
    b_escalado = b / s
        
    # Calculo el vector de traslación t resolviendo K * t = b -> t = K^-1 * b
    # Como K es triangular superior de 3x3, resuelvo por sustitución hacia atrás
    t = np.zeros(3)
    t[2] = b_escalado[2] / K[2, 2]
    t[1] = (b_escalado[1] - K[1, 2] * t[2]) / K[1, 1]
    t[0] = (b_escalado[0] - K[0, 1] * t[1] - K[0, 2] * t[2]) / K[0, 0]

    # Validación de física de la cámara: verificamos que t_z > 0 (o sea que la cámara mira hacia adelante del objeto)
    if t[2] < 0:
        # Si la traslación z queda negativa, invertimos los signos globales de la pose
        t = -t
        R = -R

    return R, t, K

# FUNCIONES EJERCICIO 3

def obtenerQ2d(imgPts, objPts):
    """
    Arma la matriz Q del sistema homogéneo Qh = 0 (ecuación 11) para estimar
    la homografía H entre puntos coplanares del mundo y puntos de la imagen.

    Parámetros:
        imgPts: matriz de n x 2 con las coordenadas (u_i, v_i) en píxeles.
        objPts: matrix de n x 2 (o n x 3) con las coordenadas (U_i, V_i) del
                plano del objeto (en mm), y W = 0.
                imgPts[i] y objPts[i] deben corresponder al mismo punto físico.

    Retorna:
        Q: array de 2n x 9, donde h = (h11, h12, h13, h21, h22, h23, h31, h32, h33)
           es el vector que debería quedar en el núcleo de Q.
    """
    imgPts = np.array(imgPts)
    objPts = np.array(objPts)

    n = len(objPts)               # cantidad de correspondencias
    Q = np.zeros((2 * n, 9))      # 2 ecuaciones por punto, 9 incógnitas

    for i in range(n):
        U, V = objPts[i, 0], objPts[i, 1]   # punto en el plano del tablero
        u, v = imgPts[i, 0], imgPts[i, 1]   # punto correspondiente en la imagen

        # Fila 2i: ecuación de u  ->  [P^T, 0^T, -u*P^T], con P = (U, V, 1)
        Q[2 * i, 0:3] = [U, V, 1]
        Q[2 * i, 6:9] = [-u * U, -u * V, -u]

        # Fila 2i+1: ecuación de v  ->  [0^T, P^T, -v*P^T]
        Q[2 * i + 1, 3:6] = [U, V, 1]
        Q[2 * i + 1, 6:9] = [-v * U, -v * V, -v]

    return Q


#----------------


def vij(H, i, j):
    """
    Arma el vector v_ij de 6 entradas tal que c_i^T B c_j = v_ij^T b,
    con b = [B11, B12, B22, B13, B23, B33]^T y c_i la i-ésima columna de H.

    Parámetros:
        H: matriz de 3 x 3 (homografía de una vista).
        i, j: índices de columna.

    Retorna:
        v: vector de 6 coordenadas.
    """
    H = np.asarray(H)
    v = np.zeros(6)
    # H[k, i] es la componente k de la columna c_i (h_{i,k+1})
    v[0] = H[0, i] * H[0, j]
    v[1] = H[0, i] * H[1, j] + H[1, i] * H[0, j]
    v[2] = H[1, i] * H[1, j]
    v[3] = H[2, i] * H[0, j] + H[0, i] * H[2, j]
    v[4] = H[2, i] * H[1, j] + H[1, i] * H[2, j]
    v[5] = H[2, i] * H[2, j]
    return v


def armarA(Hlst):
    """
    Arma la matriz A del sistema homogéneo Ab = 0, r1 y r2 ortonormales, donde
    b = [B11, B12, B22, B13, B23, B33]^T y B = K^{-T} K^{-1}.

    Parámetros:
        Hlst: lista de N matrices H_j de 3 x 3, obtenidas
              con el DLT de puntos coplanares.

    Retorna:
        A: matriz de 2N x 6. Cada vista j aporta dos filas:
           fila 2j:   v12^T          (c1^T B c2 = 0, perpendicularidad)
           fila 2j+1: (v11 - v22)^T  (c1^T B c1 = c2^T B c2, igual norma)
    """
    N = len(Hlst)                 # cantidad de vistas
    A = np.zeros((2 * N, 6))      # 2 ecuaciones por vista, 6 incógnitas en b

    for j in range(N):
        H = np.asarray(Hlst[j])
        H = H / H[2, 2]           # normalizo la escala de cada H

        v12 = vij(H, 0, 1)        # c1 con c2
        v11 = vij(H, 0, 0)        # c1 con c1
        v22 = vij(H, 1, 1)        # c2 con c2

        A[2 * j] = v12
        A[2 * j + 1] = v11 - v22

    return A


def vectorMenorSingular(M): #AUX
    """
    Devuelve el vector unitario x que minimiza ||Mx|| con ||x|| = 1.
    Es el autovector de menor autovalor de M^T M (= vector singular derecho
    del menor valor singular de M). Sirve para el DLT (Q h = 0) y para A b = 0.

    Parámetros:
        M: array de m x n.
    Retorna:
        x: array de n, de norma 1.
    """
    MtM = arroba(M.T, M)
    S, D = diagRH(MtM)          # S: autovectores (columnas), D: diagonal de autovalores
    n = MtM.shape[0]

    idx = 0
    for i in range(1, n):
        if D[i, i] < D[idx, idx]: #recorre valores singulares y se queda con el menor
            idx = i

    x = S[:, idx].copy()
    return x / norma(x, 2)


def choleskyPropio(B): #AUX
    """
    Factorización de Cholesky B = L L^T, con L triangular inferior.
    Parámetros: B (n x n) simétrica definida positiva.
    Retorna: L (n x n), o None si B no es definida positiva.
    """
    n = B.shape[0]
    L = np.zeros((n, n))
    for j in range(n):
        suma = 0.0
        for k in range(j):
            suma = suma + L[j, k] * L[j, k]
        d = B[j, j] - suma
        if d <= 0:
            return None                      # B no es definida positiva
        L[j, j] = np.sqrt(d)
        for i in range(j + 1, n):
            suma = 0.0
            for k in range(j):
                suma = suma + L[i, k] * L[j, k]
            L[i, j] = (B[i, j] - suma) / L[j, j]
    return L


def productoCruz(a, b): #AUX
    """Producto vectorial de dos vectores de R^3."""
    c = np.zeros(3)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c


def calcularMatrizK2D(A):
    """
    Calcula K a partir de A (2N x 6), donde A b = 0 con
    b = [B11, B12, B22, B13, B23, B33]^T y B = K^{-T} K^{-1}.
    Retorna K (3 x 3, triangular superior, K[2,2] = 1), o None si B no
    resulta definida positiva.
    """
    # Mínimos cuadrados homogéneos
    b = vectorMenorSingular(A)

    # Reconstruir B simétrica
    B = np.array([[b[0], b[1], b[3]],
                  [b[1], b[2], b[4]],
                  [b[3], b[4], b[5]]])

    # Signo: B debe ser definida positiva (diagonal positiva), elige signo tal que lo sea
    if B[0, 0] < 0:
        B = -1.0 * B

    # Cholesky: B = L L^T, L = K^{-T}  =>  K = (L^T)^{-1}
    L = choleskyPropio(B)
    if L is None:
        return None
    K = inversa(L.T)

    # Normalizar K[2,2] = 1
    return K / K[2, 2]


def obtenerPose2D(Hj, K):
    """
    Pose (R, t) de la vista j a partir de su homografía Hj (3x3, ya
    desnormalizada) y de K. t queda en las unidades de objPts (mm).
    """
    Hj = np.asarray(Hj, dtype=float)
    Kinv = inversa(K)

    # K^{-1} c_i  (r_i salvo escala)
    a1 = arroba(Kinv, Hj[:, 0])
    a2 = arroba(Kinv, Hj[:, 1])
    a3 = arroba(Kinv, Hj[:, 2])

    lam = 1.0 / norma(a1, 2)            # ||r1|| = 1
    r1, r2, t = lam * a1, lam * a2, lam * a3

    # Signo: el tablero debe quedar delante de la cámara (t_z > 0)
    if t[2] < 0:
        r1, r2, t = -1.0 * r1, -1.0 * r2, -1.0 * t

    # Tercer eje y matriz con columnas [r1 r2 r3]
    r3 = productoCruz(r1, r2)
    M = np.zeros((3, 3))
    M[:, 0] = r1
    M[:, 1] = r2
    M[:, 2] = r3

    # Ortonormalización con Gram-Schmidt
    R, _ = QR_con_GS(M)
    return R, t

def normalizarPuntos(pts):
    """
    Normalizacion isotropica de Hartley: lleva el centroide al origen y escala para que la
    distancia media al origen sea sqrt(2).
    Parametros: pts, array de n x 2.
    Retorna: (pts_norm de n x 2, T de 3 x 3) con p_norm = T p (en coordenadas homogeneas).
    """
    pts = np.array(pts, dtype=float)
    n = pts.shape[0]
    cx = np.mean(pts[:, 0])
    cy = np.mean(pts[:, 1])
    d = np.sqrt((pts[:, 0] - cx) ** 2 + (pts[:, 1] - cy) ** 2)
    s = np.sqrt(2.0) / np.mean(d)
    T = np.array([[s, 0.0, -s * cx],
                  [0.0, s, -s * cy],
                  [0.0, 0.0, 1.0]])
    pn = np.zeros((n, 2))
    pn[:, 0] = s * (pts[:, 0] - cx)
    pn[:, 1] = s * (pts[:, 1] - cy)
    return pn, T


def estimarHomografia(imgPts, objPts):
    Q = obtenerQ2d(imgPts, objPts)
    h = vectorMenorSingular(Q)
    H = h.reshape(3, 3)
    return H / H[2, 2]

def errorReproyeccion(R, t, K, objPts, imgPts):
    """Error medio (px) de reproyectar los puntos planos objPts (W = 0) con K [R | t]."""
    n = objPts.shape[0]
    suma = 0.0
    for i in range(n):
        P = np.array([objPts[i, 0], objPts[i, 1], 0.0])
        Xc = arroba(R, P) + t
        h = arroba(K, Xc)
        du = h[0] / h[2] - imgPts[i, 0]
        dv = h[1] / h[2] - imgPts[i, 1]
        suma = suma + np.sqrt(du * du + dv * dv)
    return suma / n

# FUNCIONES EJERCICIO 5

# FUNCIONES AUXILIARES
def extraerTiempoMs(nombreArchivo):
    """Extrae la estampa temporal en segundos con precisión de milisegundos."""
    base = os.path.basename(nombreArchivo)
    partes = base.split('.')[0].split('_')
    partes_num = [p for p in partes if p.isdigit()]
    if len(partes_num) >= 4:
        h, m, s, ms = map(float, partes_num[:4])
        return h * 3600.0 + m * 60.0 + s + (ms / 1000.0)
    return 0.0

def proyectarPuntoImagenAPiso(m, H_inv):
    """Transforma un punto en píxeles [u, v, 1] al suelo métrico [X, Y] en metros."""
    P_homog = arroba(H_inv, m)
    if abs(P_homog[2]) > 1e-15:
        return np.array([P_homog[0] / P_homog[2], P_homog[1] / P_homog[2]])
    return np.array([P_homog[0], P_homog[1]])

def estimarHomografiaHartley(imgPts, objPts):
    """
    Estima la homografía 2D H aplicando normalización isotrópica de Hartley
    y desnormalizando al final: H = T_img^{-1} * H_norm * T_obj.
    """
    img_n, T_img = normalizarPuntos(imgPts)
    obj_n, T_obj = normalizarPuntos(objPts)

    H_n = estimarHomografia(img_n, obj_n)

    T_img_inv = inversa(T_img)
    H = arroba(T_img_inv, arroba(H_n, T_obj))
    return H / H[2, 2]

# EJERCICIO 5.a: PARÁMETROS INTRÍNSECOS K (MÉTODO DE ZHANG)
def calcularKZhangRestringido(Hlst, u0=800.0, v0=600.0):
    """
    Estima K por Zhang imponiendo sesgo=0, alpha_u=alpha_v=f y centro en (u0, v0).
    Normaliza cada homografía y cada ecuación para equilibrar el ajuste lineal.
    """
    T_c = np.array([[1.0, 0.0, -u0],
                    [0.0, 1.0, -v0],
                    [0.0, 0.0, 1.0]])

    a_list, b_list = [], []
    for H in Hlst:
        Hc = arroba(T_c, H)
        Hc = Hc / norma(Hc.flatten(), 2)
        h1, h2 = Hc[:, 0], Hc[:, 1]

        a1 = h1[0] * h2[0] + h1[1] * h2[1]
        b1 = -(h1[2] * h2[2])

        a2 = (h1[0]**2 + h1[1]**2) - (h2[0]**2 + h2[1]**2)
        b2 = -(h1[2]**2 - h2[2]**2)

        n1 = np.hypot(a1, b1)
        n2 = np.hypot(a2, b2)
        a_list.extend([a1 / n1, a2 / n2])
        b_list.extend([b1 / n1, b2 / n2])

    a = np.array(a_list, dtype=float)
    b = np.array(b_list, dtype=float)

    x = np.dot(a, b) / np.dot(a, a)
    f = 1.0 / np.sqrt(x)

    return np.array([[f, 0.0, u0],
                     [0.0, f, v0],
                     [0.0, 0.0, 1.0]])


def calcularIntrinsecosSalguero(world_img_pts_csv, restringir_K=True, u0=800.0, v0=600.0):
    """
    Calcula K tratando los 5 puntos de world_img_pts.csv como coplanares (Z ~ 0)
    mediante el método de Zhang (versión libre y versión restringida).
    """
    df = pd.read_csv(world_img_pts_csv, sep=';')
    objPts_2d = df[['X', 'Y']].values

    columnas = list(df.columns)
    vistas = [c.split('_')[1] for c in columnas if c.startswith('u_') and f"v_{c.split('_')[1]}" in columnas]

    Hlst = []
    for k in range(len(vistas)):
        imgPts_k = df[[f'u_{k}', f'v_{k}']].values
        Hlst.append(estimarHomografiaHartley(imgPts_k, objPts_2d))

    K_zhang_libre = calcularMatrizK2D(armarA(Hlst))
    K_zhang_restringido = calcularKZhangRestringido(Hlst, u0=u0, v0=v0)

    if restringir_K:
        return K_zhang_restringido, K_zhang_libre
    return K_zhang_libre, K_zhang_restringido

def estimarPoseDesdeHomografia(K, H):
    """
    Extrae R (ortogonalizada en SO(3) vía SVD) y t a partir de K y la homografía H.
    """
    K_inv = inversa(K)
    K_inv_H = arroba(K_inv, H)

    h1, h2, h3 = K_inv_H[:, 0], K_inv_H[:, 1], K_inv_H[:, 2]
    n1, n2 = norma(h1, 2), norma(h2, 2)
    lam = 2.0 / (n1 + n2)

    r1, r2, t = lam * h1, lam * h2, lam * h3

    if t[2] < 0:
        r1, r2, t = -r1, -r2, -t

    r3 = productoCruz(r1, r2)
    R_raw = np.column_stack((r1, r2, r3))

    # svd_reducida retorna (U, S, V) sin trasponer
    U, _, V = svd_reducida(R_raw)
    Vt = traspuesta(V)
    R = arroba(U, Vt)

    if np.linalg.det(R) < 0:
        U[:, -1] *= -1.0
        R = arroba(U, Vt)

    C = -arroba(traspuesta(R), t)
    return R, t, C

def optimizarFocalPorReproyeccion(world_img_pts_csv, f_min=500.0, f_max=5000.0, pasos=450, u0=800.0, v0=600.0):
    """
    Busca la distancia focal f (con sesgo=0, fx=fy=f, centro=(u0, v0)) que minimiza
    el error RMS geométrico de reproyección promedio en las 4 vistas de world_img_pts.csv.
    """
    df = pd.read_csv(world_img_pts_csv, sep=';')
    objPts_2d = df[['X', 'Y']].values
    objPts_3d = df[['X', 'Y', 'Z']].values
    # Si se tratan como coplanares para la calibración:
    objPts_plano = np.column_stack((objPts_2d, np.zeros(len(objPts_2d))))

    columnas = list(df.columns)
    vistas = [c.split('_')[1] for c in columnas if c.startswith('u_') and f"v_{c.split('_')[1]}" in columnas]

    # Precalculamos las 4 homografías una sola vez
    vistas_data = []
    for k in range(len(vistas)):
        imgPts_k = df[[f'u_{k}', f'v_{k}']].values
        H_k = estimarHomografiaHartley(imgPts_k, objPts_2d)
        vistas_data.append((imgPts_k, H_k))

    mejor_f = None
    mejor_rms = float('inf')
    historial_f = np.linspace(f_min, f_max, pasos)
    historial_rms = []

    for f in historial_f:
        K_cand = np.array([[f, 0.0, u0],
                           [0.0, f, v0],
                           [0.0, 0.0, 1.0]])
        err_cuad_total = []
        for imgPts_k, H_k in vistas_data:
            R_k, t_k, _ = estimarPoseDesdeHomografia(K_cand, H_k)
            for i in range(len(objPts_plano)):
                P_cam = arroba(R_k, objPts_plano[i]) + t_k
                m_hat = arroba(K_cand, P_cam)
                u_hat, v_hat = m_hat[0] / m_hat[2], m_hat[1] / m_hat[2]
                err_cuad_total.append((imgPts_k[i, 0] - u_hat)**2 + (imgPts_k[i, 1] - v_hat)**2)

        rms = np.sqrt(np.mean(err_cuad_total))
        historial_rms.append(rms)
        if rms < mejor_rms:
            mejor_rms = rms
            mejor_f = f

    K_opt = np.array([[mejor_f, 0.0, u0],
                      [0.0, mejor_f, v0],
                      [0.0, 0.0, 1.0]])
    return K_opt, mejor_rms, historial_f, np.array(historial_rms)

# EJERCICIO 5.b: POSE DE CÁMARA (R, t) Y REPROYECCIÓN CON K[R|t]
def estimarPoseYReproyectar(file_pts, K, img_path=None, mostrar_grafico=True):
    """
    Calcula la pose (R, t) a partir de los 51 puntos sobre la calle (Z=0)
    y reproyecta directamente usando K(R * X_w + t).
    """
    data = np.loadtxt(file_pts)
    objPts = traspuesta(data[1:4, :])  # (51, 3) en metros (X, Y, Z)
    objPts[:, 2] = 0.0 # forzamos Z=0
    imgPts = traspuesta(data[4:6, :])  # (51, 2) en píxeles (u, v)
    n_pts = objPts.shape[0]

    # Homografía del plano Z=0 y pose 3D (R, t)
    H_calle = estimarHomografiaHartley(imgPts, objPts[:, :2])
    R, t, C = estimarPoseDesdeHomografia(K, H_calle)

    # 2. Reproyección directa con K[R|t]
    puntos_reproyectados = np.zeros((n_pts, 2))
    suma_err_sq = 0.0

    for i in range(n_pts):
        P_cam = arroba(R, objPts[i]) + t
        m_hat = arroba(K, P_cam)
        u_hat = m_hat[0] / m_hat[2]
        v_hat = m_hat[1] / m_hat[2]

        puntos_reproyectados[i, 0] = u_hat
        puntos_reproyectados[i, 1] = v_hat
        suma_err_sq += (imgPts[i, 0] - u_hat) ** 2 + (imgPts[i, 1] - v_hat) ** 2

    RMS = np.sqrt(suma_err_sq / float(n_pts))

    if mostrar_grafico and img_path is not None:
        im = mpimg.imread(img_path)
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.imshow(im)
        ax.scatter(imgPts[:, 0], imgPts[:, 1], c='lime', s=40, label='Puntos medidos (ground truth)', zorder=3)
        ax.scatter(puntos_reproyectados[:, 0], puntos_reproyectados[:, 1], c='red', s=35, marker='x', label='Reproyectados con K[R|t]', zorder=4)
        for i in range(n_pts):
            ax.plot([imgPts[i, 0], puntos_reproyectados[i, 0]],
                    [imgPts[i, 1], puntos_reproyectados[i, 1]], color='yellow', linewidth=0.8, zorder=2)
        ax.set_title(f'Ejercicio 5.b - Reproyección con K[R|t] (RMS: {RMS:.3f} px)')
        ax.legend()
        plt.tight_layout()
        plt.savefig('reproyeccion_salguero.png')
        plt.show()

    return R, t, C, H_calle, RMS, puntos_reproyectados

def graficarOptimizacionFocal(f_vals, rms_vals, f_opt, rms_opt, f_analitica=None):
    """Grafica el error RMS de reproyección en función de la distancia focal f."""
    plt.figure(figsize=(9, 5))
    plt.plot(f_vals, rms_vals, color='navy', linewidth=2, label='RMS geométrico (4 vistas)')
    plt.axvline(f_opt, color='red', linestyle='--', linewidth=1.5,
                label=f'f óptimo = {f_opt:.1f} px (RMS = {rms_opt:.2f} px)')
    plt.scatter([f_opt], [rms_opt], color='red', s=60, zorder=5)

    if f_analitica is not None:
        plt.axvline(f_analitica, color='darkorange', linestyle=':', linewidth=1.8,
                    label=f'f Zhang analítico = {f_analitica:.1f} px')

    plt.title('Ejercicio 5.a - Refinamiento 1D de la distancia focal f')
    plt.xlabel('Distancia focal f (px)')
    plt.ylabel('Error RMS de reproyección (px)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig('optimizacion_focal_5a.png')
    plt.show()

def graficarVelocidades(tablaVelocidades):
    """Grafica la evolución temporal de las velocidades de los tres agentes."""
    plt.figure(figsize=(10, 6))
    agentes = ['Taxi', 'Camioneta', 'Peatona']
    colores = {'Taxi': 'goldenrod', 'Camioneta': 'blue', 'Peatona': 'red'}
    estilos = {'Taxi': 'o-', 'Camioneta': 's-', 'Peatona': '^-'}

    t0 = min([f['t_i'] for f in tablaVelocidades]) if len(tablaVelocidades) > 0 else 0.0

    for agente in agentes:
        filas = [f for f in tablaVelocidades if f['agente'] == agente]
        if len(filas) > 0:
            x_vals = [(f['t_i'] + f['t_f']) / 2.0 - t0 for f in filas]
            vels = [f['v_kmh'] for f in filas]
            plt.plot(x_vals, vels, estilos[agente], color=colores[agente],
                     label=agente, linewidth=2, markersize=8)

    plt.title('Ejercicio 5.c - Evolución Temporal de Velocidades en Calle Salguero')
    plt.xlabel('Tiempo relativo desde la primera captura (s)')
    plt.ylabel('Velocidad estimada (km/h)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig('velocidades_salguero.png')
    plt.show()

def estimarVelocidadAgentes(secuencia_capturas, H_calle, tracking_agentes, factor_escala_x=1.0, factor_escala_y=1.0):
    """
    Calcula frame a frame la posición en el plano (X, Y) en metros, el desplazamiento
    y la velocidad (en m/s y km/h) del Taxi, la Camioneta y la Peatona.
    """
    H_inv = inversa(H_calle)
    Tabla_Velocidades = []
    agentes = ['Taxi', 'Camioneta', 'Peatona']

    for agente in agentes:
        if agente not in tracking_agentes:
            continue
        registros = tracking_agentes[agente]
        num_frames = min(len(registros), len(secuencia_capturas))

        for k in range(num_frames - 1):
            t_k = extraerTiempoMs(secuencia_capturas[k])
            t_k1 = extraerTiempoMs(secuencia_capturas[k + 1])
            dt = t_k1 - t_k

            u_k = registros[k][0] * factor_escala_x
            v_k = registros[k][1] * factor_escala_y
            u_k1 = registros[k + 1][0] * factor_escala_x
            v_k1 = registros[k + 1][1] * factor_escala_y

            P_k = proyectarPuntoImagenAPiso(np.array([u_k, v_k, 1.0]), H_inv)
            P_k1 = proyectarPuntoImagenAPiso(np.array([u_k1, v_k1, 1.0]), H_inv)

            dist_m = norma(P_k1 - P_k, 2)
            v_ms = dist_m / dt if dt > 0.0 else 0.0

            Tabla_Velocidades.append({
                'agente': agente,
                'frame_i': k,
                'frame_f': k + 1,
                'archivo_i': os.path.basename(secuencia_capturas[k]),
                'archivo_f': os.path.basename(secuencia_capturas[k + 1]),
                't_i': t_k,
                't_f': t_k1,
                'X_i': P_k[0],
                'Y_i': P_k[1],
                'X_f': P_k1[0],
                'Y_f': P_k1[1],
                'dist_m': dist_m,
                'dt': dt,
                'v_ms': v_ms,
                'v_kmh': v_ms * 3.6
            })

    return Tabla_Velocidades

# Tracking manual (u, v) medido en imágenes de 1600x1200, ordenado cronológicamente
tracking_agentes = {
        'Peatona': [
            (420.0,500.0),  # 10:18:47.351 - cruzando senda peatonal
            (434.0,502.0),  # 10:18:47.461
            (456.0, 500.0),  # 10:18:47:554
            (470.0, 500.0),  # 10:18:47.679
            (478.0, 500.0),  # 10:18:47.835
            (498.0, 500.0),  # 10:18:47.882
            (514.0, 499.0),  # 10:18:48.022
            (521.0, 497.0),  # 10:18:48.163
            (535.0, 496.0),  # 10:18:48.209
            (546.0, 496.0),  # 10:18:48.350
            (554.0, 498.0),  # 10:18:48.475
            (584.0, 495.0),  # 10:18:48.724
            (597.0, 497.0),  # 10:18:48.818
            (619.0, 499.0),  # 10:18:49.005
            (640.0, 497.0),  # 10:18:49.099
            (647.0, 495.0),  # 10:18:49.286
            (668.0, 494.0),  # 10:18:49.317
            (686.0, 497.0),  # 10:18:49.457
            (694.0, 495.0),  # 10:18:49.582
            (700.0, 495.0),  # 10:18:49.660
            (728.0, 494.0),  # 10:18:49.785
            (741.0, 493.0),  # 10:18:49.879
            (743.0, 495.0),  # 10:18:49.988
            (760.0, 495.0),  # 10:18:50.128
            (772.0, 495.0),  # 10:18:50.206
            (785.0, 494.0),  # 10:18:50.300
            (792.0, 495.0),  # 10:18:50.471
            (808.0, 494.0),  # 10:18:50.534
            (821.0, 494.0),  # 10:18:50.674
            (834.0, 494.0),  # 10:18:50.783
            (844.0, 491.0),  # 10:18:50.893
            (852.0, 496.0),  # 10:18:50.986
            (863.0, 494.0),  # 10:18:51.173
            (877.0, 494.0),  # 10:18:51.251
            (888.0, 492.0),  # 10:18:51.345
            (904.0, 491.0),  # 10:18:51.470
            (909.0, 492.0),  # 10:18:51.563
            (918.0, 492.0),  # 10:18:51.657
            (930.0, 491.0),  # 10:18:51.844
            (944.0, 490.0),  # 10:18:51.875
            (958.0, 492.0),  # 10:18:52.000
            (964.0, 492.0),  # 10:18:52.109
            (973.0, 492.0),  # 10:18:52.219
            (978.0, 492.0),  # 10:18:52.312
            (980.0, 492.0),  # 10:18:52.453
            (992.0, 492.0),  # 10:18:52.546
        ],
        'Taxi': [
            (930.0, 269.0),  # 10:18:47.351 - carril izquierdo acercándose
            (929.0, 271.0),  # 10:18:47.461
            (923.0, 272.0),  # 10:18:47:554
            (920.0, 275.0),  # 10:18:47.679
            (915.0, 275.0),  # 10:18:47.835
            (912.0, 277.0),  # 10:18:47.882
            (907.0, 281.0),  # 10:18:48.022
            (905.0, 280.0),  # 10:18:48.163
            (900.0, 284.0),  # 10:18:48.209
            (897.0, 285.0),  # 10:18:48.350
            (893.0, 286.0),  # 10:18:48.475
            (886.0, 290.0),  # 10:18:48.724
            (879.0, 292.0),  # 10:18:48.818
            (867.0, 296.0),  # 10:18:49.005
            (861.0, 299.0),  # 10:18:49.099
            (857.0, 301.0),  # 10:18:49.286
            (844.0, 303.0),  # 10:18:49.317
            (840.0, 308.0),  # 10:18:49.457
            (833.0, 310.0),  # 10:18:49.582
            (829.0, 312.0),  # 10:18:49.660
            (816.0, 317.0),  # 10:18:49.785
            (806.0, 324.0),  # 10:18:49.879
            (802.0, 324.0),  # 10:18:49.988
            (792.0, 331.0),  # 10:18:50.128
            (780.0, 335.0),  # 10:18:50.206
            (770.0, 340.0),  # 10:18:50.300
            (762.0, 345.0),  # 10:18:50.471
            (751.0, 352.0),  # 10:18:50.534
            (737.0, 357.0),  # 10:18:50.674
            (729.0, 361.0),  # 10:18:50.783
            (714.0, 366.0),  # 10:18:50.893
            (704.0, 370.0),  # 10:18:50.986
            (684.0, 384.0),  # 10:18:51.173
            (663.0, 390.0),  # 10:18:51.251
            (639.0, 404.0),  # 10:18:51.345
            (630.0, 408.0),  # 10:18:51.470
            (604.0, 418.0),  # 10:18:51.563
            (589.0, 425.0),  # 10:18:51.657
            (561.0, 442.0),  # 10:18:51.844
            (521.0, 456.0),  # 10:18:51.875
            (486.0, 477.0),  # 10:18:52.000
            (463.0, 484.0),  # 10:18:52.109
            (416.0, 503.0),  # 10:18:52.219
            (384.0, 516.0),  # 10:18:52.312
            (323.0, 547.0),  # 10:18:52.453
            (284.0, 557.0),  # 10:18:52.546
        ],
        'Camioneta': [
            (995.0, 261.0),  # 10:18:47.351 - carril central detrás del taxi
            (994.0, 260.0),  # 10:18:47.461
            (994.0, 263.0),  # 10:18:47:554
            (993.0, 263.0),  # 10:18:47.679
            (991.0, 266.0),  # 10:18:47.835
            (990.0, 266.0),  # 10:18:47.882
            (991.0, 268.0),  # 10:18:48.022
            (988.0, 269.0),  # 10:18:48.163
            (989.0, 271.0),  # 10:18:48.209
            (990.0, 270.0),  # 10:18:48.350
            (989.0, 272.0),  # 10:18:48.475 
            (986.0, 275.0),  # 10:18:48.724
            (983.0, 278.0),  # 10:18:48.818
            (983.0, 281.0),  # 10:18:49.005
            (981.0, 281.0),  # 10:18:49.099
            (979.0, 285.0),  # 10:18:49.286
            (975.0, 285.0),  # 10:18:49.317
            (972.0, 287.0),  # 10:18:49.457
            (972.0, 289.0),  # 10:18:49.582
            (971.0, 291.0),  # 10:18:49.660
            (964.0, 296.0),  # 10:18:49.785
            (962.0, 298.0),  # 10:18:49.879
            (960.0, 300.0),  # 10:18:49.988
            (957.0, 303.0),  # 10:18:50.128
            (952.0, 307.0),  # 10:18:50.206
            (948.0, 311.0),  # 10:18:50.300
            (945.0, 311.0),  # 10:18:50.471
            (941.0, 316.0),  # 10:18:50.534
            (934.0, 321.0),  # 10:18:50.674
            (932.0, 324.0),  # 10:18:50.783
            (928.0, 328.0),  # 10:18:50.893
            (925.0, 332.0),  # 10:18:50.986
            (915.0, 336.0),  # 10:18:51.173
            (912.0, 344.0),  # 10:18:51.251
            (900.0, 350.0),  # 10:18:51.345
            (896.0, 355.0),  # 10:18:51.470
            (886.0, 363.0),  # 10:18:51.563
            (882.0, 368.0),  # 10:18:51.657
            (873.0, 374.0),  # 10:18:51.844
            (862.0, 386.0),  # 10:18:51.875
            (848.0, 398.0),  # 10:18:52.000
            (841.0, 403.0),  # 10:18:52.109
            (824.0, 415.0),  # 10:18:52.219
            (818.0, 421.0),  # 10:18:52.312
            (795.0, 438.0),  # 10:18:52.453
            (788.0, 449.0),  # 10:18:52.546
        ]
    }

def suavizarTrackingManual(tracking_dict, secuencia_capturas, H_calle, grado=2):
    """
    Suaviza el ruido de clic manual (+/- 1-2 px) ajustando la trayectoria en el plano
    del suelo (X(t), Y(t)) en función del tiempo REAL de cada captura mediante
    regresión polinomial (QR_con_GS y res_tri de py).
    """
    H_inv_temp = inversa(H_calle)
    tracking_suave = {}

    # Tiempos reales de las capturas
    tiempos = np.array([extraerTiempoMs(f) for f in secuencia_capturas], dtype=float)

    for ag, pts in tracking_dict.items():
        n = min(len(pts), len(tiempos))
        t_ag = tiempos[:n]
        # Normalizamos el tiempo a [0, 1] para buen condicionamiento de Vandermonde
        t_norm = (t_ag - t_ag[0]) / (t_ag[-1] - t_ag[0]) if (t_ag[-1] - t_ag[0]) > 0 else np.linspace(0.0, 1.0, n)

        # Proyectamos los clics al piso (Z=0) en metros
        XY = np.zeros((n, 2))
        for k in range(n):
            XY[k] = proyectarPuntoImagenAPiso(np.array([pts[k][0], pts[k][1], 1.0]), H_inv_temp)

        # Matriz de Vandermonde V (n x (grado + 1)) en función de t_norm
        V = np.zeros((n, grado + 1))
        for d in range(grado + 1):
            V[:, d] = t_norm ** d

        Q_v, R_v = QR_con_GS(V)
        coef_X = res_tri(R_v, arroba(traspuesta(Q_v), XY[:, 0]), inferior=False)
        coef_Y = res_tri(R_v, arroba(traspuesta(Q_v), XY[:, 1]), inferior=False)

        X_suave = arroba(V, coef_X)
        Y_suave = arroba(V, coef_Y)

        # Reproyectamos a coordenadas de imagen (u, v) con H_calle
        pts_suaves = []
        for k in range(n):
            m_hat = arroba(H_calle, np.array([X_suave[k], Y_suave[k], 1.0]))
            pts_suaves.append((float(m_hat[0] / m_hat[2]), float(m_hat[1] / m_hat[2])))

        tracking_suave[ag] = pts_suaves

    return tracking_suave