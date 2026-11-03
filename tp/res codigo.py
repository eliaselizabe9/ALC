import numpy as np
from alc import arroba, diagRH, QR_con_HH

# EJERCICIO 2.a
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


# EJERCICIO 2.b
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
    return M


# EJERCICIOS 2.c y 2.d: DESCOMPOSICIÓN DE M EN K, R Y t
def descomposicion_RQ_3x3(A):
    """
    Realizo la descomposición RQ de una matriz A de 3x3 mediante la factorización QR de alc.
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
    
    # Normalizo K para asegurar que el elemento homogéneo K[2,2] valga 1
    K = K / K[2, 2]
    
    # Impongo las restricciones opcionales sobre la matriz de intrínsecos K
    if fijar_sesgo_cero:
        K[0, 1] = 0.0  # gamma = 0
        
    if igualar_focales:
        alpha_prom = (K[0, 0] + K[1, 1]) / 2.0
        K[0, 0] = alpha_prom
        K[1, 1] = alpha_prom
        
    # Calculo el vector de traslación t resolviendo K * t = b -> t = K^-1 * b
    # Como K es triangular superior de 3x3, resuelvo por sustitución hacia atrás
    t = np.zeros(3)
    t[2] = b[2] / K[2, 2]
    t[1] = (b[1] - K[1, 2] * t[2]) / K[1, 1]
    t[0] = (b[0] - K[0, 1] * t[1] - K[0, 2] * t[2]) / K[0, 0]
    
    return R, t, K