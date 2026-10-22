import numpy as np

# 1. NORMA VECTORIAL
def norma(x, p):
    """Calcula la p-norma de un vector x."""

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


# 2. NORMALIZACIÓN DE VECTORES
def normaliza(X, p):
    """Devuelve una lista con los vectores de X normalizados."""

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


# 3. ESTIMACIÓN DE NORMA MATRICIAL INDUCIDA (MONTE CARLO)
def normaMatMC(A, q, p, Np):
    """Estima la norma matricial inducida y devuelve (máximo, x_max)."""

    if A.ndim != 2 or A.size == 0 or Np <= 0:
        return None

    m, n = A.shape
    max_norma = -1.0
    x_max = np.zeros(n, dtype=float)

    for _ in range(Np):
        v_rand = np.random.normal(0.0, 1.0, size=n)

        vectores_norm = normaliza([v_rand], p)
        if vectores_norm is None:
            return None

        x_unit = vectores_norm[0]

        y = np.zeros(m, dtype=float)

        for i in range(m):
            suma_fila = 0.0
            for j in range(n):
                suma_fila += A[i, j] * x_unit[j]
            y[i] = suma_fila

        val_q = norma(y, q)

        if val_q is None:
            return None

        if val_q > max_norma:
            max_norma = val_q
            x_max = x_unit.copy()

    return (float(max_norma), x_max)


# 4. NORMA MATRICIAL EXACTA (1 O INFINITO)
def normaExacta(A, p=[1,'inf']):
    """Calcula la norma matricial subordinada para p=1 o p='inf'."""

    if A.ndim != 2 or A.size == 0:
        return None

    if p == 1:
        return norma1(A)

    if p == 'inf':
        return normaInf(A)

    return None


def norma1(A):
    """Máxima suma de valores absolutos por columna."""

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
    """Máxima suma de valores absolutos por fila."""

    m, n = A.shape
    max_fila = 0.0

    for i in range(m):
        suma_fila = 0.0
        for j in range(n):
            suma_fila += abs(A[i, j])

        if suma_fila > max_fila:
            max_fila = suma_fila

    return float(max_fila)


# 5. NÚMERO DE CONDICIÓN ESTIMADO (MONTE CARLO)
def condMC(A, p, Np=1000):
    """Estima el número de condición mediante Monte Carlo."""

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


# 6. INVERSA POR GAUSS-JORDAN CON PIVOTEO PARCIAL
def inversa(A):
    """Calcula A^(-1) mediante Gauss-Jordan con pivoteo parcial."""

    if A.ndim != 2 or A.shape[0] != A.shape[1] or A.size == 0:
        return None

    n = A.shape[0]
    aug = np.hstack([A.copy(), np.eye(n, dtype=float)])

    for col in range(n):
        max_fila = col + int(np.argmax(np.abs(aug[col:, col])))
        max_val = abs(aug[max_fila, col])

        if max_val < 1e-14:
            return None

        if max_fila != col:
            aug[[col, max_fila]] = aug[[max_fila, col]]

        pivote = aug[col, col]
        aug[col, :] /= pivote

        for fila in range(n):
            if fila != col:
                factor = aug[fila, col]
                aug[fila, :] -= factor * aug[col, :]

    return aug[:, n:].copy()


# 7. NÚMERO DE CONDICIÓN EXACTO (PARA p=1 O INFINITO)
def condExacta(A, p):
    """Calcula Cond_p(A) mediante normas matriciales exactas."""

    A_inv = inversa(A)

    if A_inv is None:
        return None

    norma_A = normaExacta(A, p)
    norma_Ainv = normaExacta(A_inv, p)

    if norma_A is None or norma_Ainv is None:
        return None

    return float(norma_A * norma_Ainv)
