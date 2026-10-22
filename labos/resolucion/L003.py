import numpy as np

# 1. NORMA VECTORIAL
def norma(x, p):
    """Calcula la p-norma de un vector x en R^n utilizando NumPy[cite: 1].
    
    Parámetros:
        x (array_like): Vector de entrada[cite: 1].
        p (int, float, str): Parámetro en [1, inf) o 'inf'/'Inf'[cite: 1].
    """
    x = np.asarray(x, dtype=float).ravel()
    n = len(x)
    if n == 0:
        return None

    if isinstance(p, str) and p.lower() == "inf":
        val_max = 0.0
        for i in range(n):
            val_abs = abs(x[i])
            if val_abs > val_max:
                val_max = val_abs
        return float(val_max)
    else:
        if p < 1:
            return None

        suma = 0.0
        for i in range(n):
            suma += abs(x[i]) ** p

        return float(suma ** (1.0 / p))


# 2. NORMALIZACIÓN DE VECTORES
def normaliza(X, p):
    """Devuelve una lista con los vectores de X normalizados (norma p igual a 1)[cite: 1].
    
    Parámetros:
        X (list de array_like): Conjunto de vectores[cite: 1].
        p: Parámetro de norma[cite: 1].
    """
    X_norm = []

    for k, v in enumerate(X):
        v_arr = np.asarray(v, dtype=float).ravel()
        norma_v = norma(v_arr, p)

        if norma_v <= 1e-15:
            return None

        # Creación y escalado de las componentes
        n = len(v_arr)
        v_unitario = np.zeros(n, dtype=float)
        for i in range(n):
            v_unitario[i] = v_arr[i] / norma_v

        X_norm.append(v_unitario)

    return X_norm


# 3. ESTIMACIÓN DE NORMA MATRICIAL INDUCIDA (MONTE CARLO)
def normaMatMC(A, q, p, Np):
    """Estima ||A||_{q, p} y el vector unitario x_max mediante Monte Carlo[cite: 1].
    
    Restricción: El producto matriz-vector se calcula componente a componente[cite: 1].
    """
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.size == 0:
        return None

    m, n = A.shape
    max_norma = -1.0
    x_max = np.zeros(n, dtype=float)

    for _ in range(Np):
        # Muestreo estocástico ~ Normal(0, 1)[cite: 1]
        v_rand = np.random.normal(0.0, 1.0, size=n)

        # Proyección sobre la esfera unitaria de norma p[cite: 1]
        x_unit = normaliza([v_rand], p)[0]

        # Multiplicación matriz-vector explícita sin operadores directos (@ / np.matmul)[cite: 1]
        y = np.zeros(m, dtype=float)
        for i in range(m):
            suma_fila = 0.0
            for j in range(n):
                suma_fila += A[i, j] * x_unit[j]
            y[i] = suma_fila

        # Evaluación de la q-norma[cite: 1]
        val_q = norma(y, q)

        if val_q > max_norma:
            max_norma = val_q
            x_max = x_unit.copy()

    return (float(max_norma), x_max)


# 4. NORMA MATRICIAL EXACTA (1 O INFINITO)
def normaExacta(A, p=[1,'inf']):
    """Calcula analíticamente la norma matricial subordinada para p=1 o p='inf'[cite: 1].
    
    Parámetros:
        A (array_like): Matriz bidimensional en R^{m x n}[cite: 1].
        p (int, str): 1 para suma máxima por columnas, 'inf'/'Inf' para suma máxima por filas[cite: 1].
        
    Retorna:
        float: Valor de la norma calculada[cite: 1].
    """
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.size == 0:
        return None

    m, n = A.shape

    if p == 1:
        return norma1(A)

    elif p == "inf":
        return normaInf(A)

    else:
        return None

# FUNCIONES AUXILIARES: norma1, normaInf
def norma1(A):
    max_col = 0.0
    for j in range(n):
        suma_col = 0.0
        for i in range(m):
            suma_col += abs(A[i, j])
        if suma_col > max_col:
            max_col = suma_col
    return float(max_col)

def normaInf(A):
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
    """Estima Cond_p(A) utilizando Monte Carlo para calcular ||A||_p y ||A^{-1}||_p[cite: 1]."""
    A_inv = inversa(A)
    if A_inv is None:
        return None

    res_A = normaMatMC(A, p, p, Np)
    norma_A = res_A[0]  # Elemento max_norma[cite: 1]

    res_Ainv = normaMatMC(A_inv, p, p, Np)
    norma_Ainv = res_Ainv[0]

    return float(norma_A * norma_Ainv)

# FUNCIÓN AUXILIAR: inversa
def inversa(A):
    """Calcula A^{-1} mediante reducción por Gauss-Jordan con pivoteo parcial."""
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None

    n = A.shape[0]
    # Matriz aumentada [A | I]
    aug = np.hstack([A.copy(), np.eye(n, dtype=float)])

    for col in range(n):
        # Pivoteo parcial
        max_fila = col + int(np.argmax(np.abs(aug[col:, col])))
        max_val = abs(aug[max_fila, col])

        if max_val < 1e-14:
            return None  # Singular o no invertible[cite: 1]

        # Intercambio de filas
        if max_fila != col:
            aug[[col, max_fila]] = aug[[max_fila, col]]

        # Normalización de la fila pivote
        pivote = aug[col, col]
        aug[col, :] /= pivote

        # Eliminación del resto de las filas
        for fila in range(n):
            if fila != col:
                factor = aug[fila, col]
                aug[fila, :] -= factor * aug[col, :]

    return aug[:, n:].copy()


# 6. NÚMERO DE CONDICIÓN EXACTO
def condExacta(A, p):
    """Calcula analíticamente Cond_p(A) para p en {1, 'inf'}[cite: 1]."""
    A_inv = inversa(A)
    if A_inv is None:
        return None

    # Al retornar directamente un flotante, se multiplican de forma directa[cite: 1]
    norma_A = normaExacta(A, p)
    norma_Ainv = normaExacta(A_inv, p)

    return float(norma_A * norma_Ainv)