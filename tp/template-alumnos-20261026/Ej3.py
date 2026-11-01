''' Ejercicio 3: Resolver la calibracion para el caso de puntos copalanares. Las imagenes muestran un tablero de ajedrez
 plano visto desde distintas poses, con 12x8 cuadrados, es decir 11x7 esquinas internas, y cada cyadrado de 30x30mm.
 Tomar un sistema de coordenadas del tablero un origen en una de las esquinas internas, con U y V a lo largo y W = 0,
 de modo que las esquinas internas coordenadas (30i,30j)mm, con i ={0,...,10} y j = {0,...,6}.'''


import numpy as np

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

