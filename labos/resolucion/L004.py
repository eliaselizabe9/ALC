import numpy as np

# 1. DESCOMPOSICIÓN LU
def calculaLU(A):
    """
    Calcula la factorizacion LU sin pivoteo tal que A = LU.
    L es triangular inferior con 1s en la diagonal y U es triangular superior.
    Retorna (L, U, ops) o None si la matriz no admite la factorizacion.
    """
    A_arr = np.array(A, dtype=float)
    if A_arr.ndim != 2:
        return None

    n, m = A_arr.shape
    if n != m or n == 0:
        return None

    L = np.eye(n, dtype=float)
    U = A_arr.copy()
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

# 2. RESOLVER TRIANGULAR
def res_tri(T, b, inferior=True):
    """
    Resuelve el sistema lineal Tx = b donde T es una matriz triangular.
    Si inferior=True, realiza sustitucion hacia adelante.
    Si inferior=False, realiza sustitucion hacia atras.
    """
    T_arr = np.array(T, dtype=float)
    b_arr = np.array(b, dtype=float).flatten()
    n = len(b_arr)

    x = np.zeros(n, dtype=float)

    if inferior:
        for i in range(n):
            if abs(T_arr[i, i]) <= 1e-15:
                return None

            suma = 0.0
            for j in range(i):
                suma += T_arr[i, j] * x[j]

            x[i] = (b_arr[i] - suma) / T_arr[i, i]
    else:
        for i in range(n - 1, -1, -1):
            if abs(T_arr[i, i]) <= 1e-15:
                return None

            suma = 0.0
            for j in range(i + 1, n):
                suma += T_arr[i, j] * x[j]

            x[i] = (b_arr[i] - suma) / T_arr[i, i]

    return x


# 3. INVERSA
def inversa(A):
    """
    Calcula la inversa de A resolviendo n sistemas triangulares dobles a partir de LU.
    Retorna la matriz inversa X o None si no es invertible / no admite LU.
    """
    res = calculaLU(A)
    if res is None:
        return None

    L, U, _ = res
    A_arr = np.array(A, dtype=float)
    n = A_arr.shape[0]
    X = np.zeros((n, n), dtype=float)

    for j in range(n):
        e_j = np.zeros(n, dtype=float)
        e_j[j] = 1.0

        y_j = res_tri(L, e_j, inferior=True)
        x_j = res_tri(U, y_j, inferior=False)

        for i in range(n):
            X[i, j] = x_j[i]

    return X


# 4. CALCULAR LDV
def calculaLDV(A):
    """
    Calcula la factorizacion LDV tal que A = L * D * V.
    L es triangular inferior unitaria, D es diagonal y V es triangular superior unitaria.
    Retorna (L, D, V) o None si falla.
    """
    res = calculaLU(A)
    if res is None:
        return None

    L, U, _ = res
    A_arr = np.array(A, dtype=float)
    n = A_arr.shape[0]

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


# 5. SIMÉTRICA DEFINIDA POSITIVA
def esSDP(A, atol=1e-8):
    """
    Determina si una matriz A es Simetrica Definida Positiva (SDP).
    Verifica simetria bajo tolerancia atol y evalua positividad estricta
    de los elementos de la diagonal D mediante calculaLDV.
    """
    A_arr = np.array(A, dtype=float)
    if A_arr.ndim != 2:
        return False

    n, m = A_arr.shape
    if n != m or n == 0:
        return False

    # Verificacion explicita de simetria
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A_arr[i, j] - A_arr[j, i]) > atol:
                return False

    # Obtencion de factorizacion LDV
    res = calculaLDV(A_arr)
    if res is None:
        return False

    L, D, V = res

    # Verificacion de elementos diagonales positivos
    for i in range(n):
        if D[i, i] <= atol:
            return False

    return True