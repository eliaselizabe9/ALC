import numpy as np

# 1. Fact. QR con Gram-Schmidt
def QR_con_GS(A, tol=1e-12, retorna_nops=False):
    """
    A una matriz de n x n 
    tol la tolerancia con la que se filtran elementos nulos en R
    retorna_nops permite (opcionalmente) retornar el numero de operaciones realizado
    retorna matrices Q y R calculadas con Gram Schmidt (y como tercer argumento opcional, el numero de operaciones).
    Si la matriz A no es de n x n, debe retornar None
    """
    if A.ndim != 2:
        return None
    
    n, m = A.shape
    if n != m:
        return None

    V = A.copy()
    Q = np.zeros((n, n), dtype=float)
    R = np.zeros((n, n), dtype=float)
    ops = 0

    for i in range(n):
        # Calculamos la norma euclídea de la columna V[:, i]
        norma_sq = 0.0
        for k in range(n):
            norma_sq += V[k, i] * V[k, i]
            ops += 2

        r_ii = norma_sq**(1/2)
        ops += 1 

        if r_ii < tol:
            R[i, i] = 0.0
            for k in range(n):
                Q[k, i] = 0.0
        else:
            R[i, i] = r_ii
            for k in range(n):
                Q[k, i] = V[k, i] / r_ii
                ops += 1

        # Ortogonalización sobre las columnas restantes (j = i+1 ... n-1)
        for j in range(i + 1, n):
            dot = 0.0
            for k in range(n):
                dot += Q[k, i] * V[k, j]
                ops += 2  # 1 mult + 1 suma

            if abs(dot) < tol:
                dot = 0.0

            R[i, j] = dot

            for k in range(n):
                V[k, j] -= dot * Q[k, i]
                ops += 2 

    if retorna_nops:
        return Q, R, ops
    return Q, R


# 2. Fact. QR con Houseolder
def QR_con_HH(A, tol=1e-12, extras=False):
    """
    A una matriz de m x n (m>=n)
    tol la tolerancia con la que se filtran elementos nulos en R
    retorna matrices Q y R calculadas con reflexiones de Householder
    Si la matriz A no cumple m>=n, debe retornar None
    extras : bool, opcional
        Si es True, devuelve informacion extra sobre el proceso de factorizacion.
        Por defecto es False. Esto lo hacemos para poder graficar el proceso.
    Devuelve la factorizacion QR de A usando reflectores de Householder.
    Devuelve: 
        Q, R, extra_info (si extras es True)
        Q, R (si extras es False)
    extra_info es un diccionario con la clave:
        'R_matrices': lista de las matrices R en cada paso
        'Q_matrices': lista de las matrices Q en cada paso
    """
    if A.ndim != 2:
        return None

    m, n = A.shape
    if m < n:
        return None

    R = A.copy()
    Q = np.eye(m, dtype=float)

    R_lista = []
    Q_lista = []

    # Se itera hasta min(m - 1, n) para matrices cuadradas o rectangulares
    num_pasos = min(m - 1, n) if m == n else n

    for k in range(num_pasos):
        sub_len = m - k

        norma_x_sq = 0.0
        for i in range(k, m):
            norma_x_sq += R[i, k] * R[i, k]
        norma_x = np.sqrt(norma_x_sq)

        if norma_x > tol:
            # Signo de R[k, k] (si es 0.0 se toma 1.0)
            signo = 1.0 if R[k, k] >= 0.0 else -1.0
            alpha = -signo * norma_x

            # Construcción de v
            v = np.zeros(sub_len, dtype=float)
            v[0] = R[k, k] - alpha
            for i in range(1, sub_len):
                v[i] = R[k + i, k]

            norma_v_sq = 0.0
            for i in range(sub_len):
                norma_v_sq += v[i] * v[i]
            norma_v = np.sqrt(norma_v_sq)

            if norma_v > 1e-15:
                # Normalización u = v / ||v||
                u = np.zeros(sub_len, dtype=float)
                for i in range(sub_len):
                    u[i] = v[i] / norma_v

                # Actualización R = H_k * R (submatriz R[k:m, k:n])
                for j in range(k, n):
                    dot = 0.0
                    for i in range(sub_len):
                        dot += u[i] * R[k + i, j]

                    for i in range(sub_len):
                        R[k + i, j] -= 2.0 * u[i] * dot

                # Forzar el elemento diagonal exacto y ceros teóricos debajo en la columna k
                R[k, k] = alpha
                for i in range(k + 1, m):
                    R[i, k] = 0.0

                # Actualización Q = Q * H_k (subcolumnas Q[:, k:m])
                for i in range(m):
                    dot = 0.0
                    for j in range(sub_len):
                        dot += Q[i, k + j] * u[j]

                    for j in range(sub_len):
                        Q[i, k + j] -= 2.0 * dot * u[j]

        if extras:
            R_lista.append(R.copy())
            Q_lista.append(Q.copy())

    # Filtrado final de tolerancia
    for i in range(m):
        for j in range(n):
            if i > j or abs(R[i, j]) < tol:
                R[i, j] = 0.0

    if extras:
        extra_info = {'R_matrices': R_lista, 'Q_matrices': Q_lista}
        return Q, R, extra_info
    return Q, R


# 3. Fact. QR
def calculaQR(A, metodo='RH', tol=1e-12):
    """
    A una matriz de n x n 
    tol la tolerancia con la que se filtran elementos nulos en R    
    metodo = ['RH','GS'] usa reflectores de Householder (RH) o Gram Schmidt (GS) para realizar la factorizacion
    retorna matrices Q y R calculadas con Gram Schmidt (y como tercer argumento opcional, el numero de operaciones)
    Si el metodo no esta entre las opciones, retorna None
    """
    if A.ndim != 2:
        return None

    n, m = A.shape
    if n != m:
        return None

    if metodo == 'RH':
        return QR_con_HH(A, tol=tol, extras=False)
    elif metodo == 'GS':
        return QR_con_GS(A, tol=tol, retorna_nops=False)
    else:
        return None
