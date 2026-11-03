''' Ejercicio 3: Resolver la calibracion para el caso de puntos copalanares. Las imagenes muestran un tablero de ajedrez
 plano visto desde distintas poses, con 12x8 cuadrados, es decir 11x7 esquinas internas, y cada cyadrado de 30x30mm.
 Tomar un sistema de coordenadas del tablero un origen en una de las esquinas internas, con U y V a lo largo y W = 0,
 de modo que las esquinas internas coordenadas (30i,30j)mm, con i ={0,...,10} y j = {0,...,6}.'''


import numpy as np
import alc

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
    MtM = alc.arroba(M.T, M)
    S, D = alc.diagRH(MtM)          # S: autovectores (columnas), D: diagonal de autovalores
    n = MtM.shape[0]

    idx = 0
    for i in range(1, n):
        if D[i, i] < D[idx, idx]: #recorre valores singulares y se queda con el menor
            idx = i

    x = S[:, idx].copy()
    return x / alc.norma(x, 2)


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
    K = alc.inversa(L.T)

    # Normalizar K[2,2] = 1
    return K / K[2, 2]


def obtenerPose2D(Hj, K):
    """
    Pose (R, t) de la vista j a partir de su homografía Hj (3x3, ya
    desnormalizada) y de K. t queda en las unidades de objPts (mm).
    """
    Hj = np.asarray(Hj, dtype=float)
    Kinv = alc.inversa(K)

    # K^{-1} c_i  (r_i salvo escala)
    a1 = alc.arroba(Kinv, Hj[:, 0])
    a2 = alc.arroba(Kinv, Hj[:, 1])
    a3 = alc.arroba(Kinv, Hj[:, 2])

    lam = 1.0 / alc.norma(a1, 2)            # ||r1|| = 1
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
    R, _ = alc.QR_con_GS(M)
    return R, t