'''Laboratorio 6: Metodo de la potencia'''

import numpy as np

def norma2(v): #calcula la norma 2 de un vector
    v = np.array(v)
    n = len(v)
    contador = 0

    for i in range(n):
        contador+= v[i]**2

    return contador**(1/2)

def f(A, v):
    """f_A(v) = Av / ||Av||"""
    w = A @ v
    return w / norma2(w)

def f2(A, v):
    """f_A^2(v) = f_A(f_A(v))"""
    return f(A, f(A, v))

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


def householder(v1):
    n = len(v1)
    e1 = np.zeros(n)
    e1[0] = 1.0
    u = e1 - v1
    nu2 = u @ u
    if nu2 < 1e-30:                  # v1 ya es e1: H = I (evita dividir por 0)
        return np.eye(n)
    return np.eye(n) - 2.0 * np.outer(u, u) / nu2

def diagRH(A, tol=1e-15, K=1000):
    '''A una matriz de nxn.
    tol la tolerancia en la diferencia entre un paso y el siguiente de la estimacion del autovector.
    K el numero maximo de iteraciones a realizarse.
    Retorna una matriz de autovectores S y matrix de autovalores D, tal que A = SDS_t.
    Si la matriz A no es simetrica debe retornar None'''
    n = A.shape[0]
    if n == 1:                       # caso base 
        return np.eye(1), A.copy()

    v1, lam1, _ = metpot2k(A, tol, K)
    H = householder(v1)
    B = H @ A @ H.T

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
    S = H @ S_1
    return S, D

# -------------AUX----------------

def norma1(A): #aux
    A = np.array(A)
    n,m = np.shape(A)
    normas = np.zeros(m)  
    for j in range(m):
        suma_columna = 0
        for i in range(n):
            suma_columna += abs(A[i, j])
        normas[j] = suma_columna
            
    return np.max(normas)

def normaInf(A):  #aux
    A = np.array(A)
    n,m = np.shape(A)
    normas = np.zeros(n)  
    for i in range(n):
        suma_fila = 0
        for j in range(m):
            suma_fila += abs(A[i, j])
        normas[i] = suma_fila
            
    return np.max(normas)

# Norma exacta
def normaExacta(A, p=[1,'inf']):
    if p == 1:
        return norma1(A)
    elif p == 'inf':
        return normaInf(A)
    elif p == [1,'inf']:
        return [norma1(A), normaInf(A)]
    else:
        return None




    
