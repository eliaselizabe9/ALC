# Labo 0 



#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

# Labo 1

def error(x,y):
    '''
    Recibe dos numeros x e y, y calcula el error de aproximar x usando y en float64
    '''
    x_f64 = np.float64(x)
    y_f64 = np.float64(y)
    
    err = abs(x_f64 - y_f64)
    
    return err


def error relativo(x,y):
    '''
    Recibe dos numeros x e y, y calcula el error relativo de aproximar x usando y en float64
    '''
    x_f64 = np.float64(x)
    y_f64 = np.float64(y)
    
    if x_f64 == 0:
        raise ValueError("El valor de x no puede ser cero para calcular el error relativo.")
        
    return abs(x_f64 - y_f64) / abs(x_f64)


def matricesIguales(A,B):
    '''
    Devuelve True si ambas matrices son iguales y False en otro caso.
    Considerar que las matrices pueden tener distintas dimensiones, ademas de distintos valores.
    '''
    atol=10**(-12)
    rtol=10**(-9)
    if A.shape != B.shape:
        return False
    if A.shape == B.shape:
        columnas=range(A.shape[0])
        filas= range(A.shape[1])
        for fil in filas:
            for col in columnas:
                b=float(B[fil][col]) 
                a= float(A[fil][col])
                if error(a, b) >= atol + rtol*(abs(b)):
                    return False
    return True

#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

# Labo 2

def rota(theta):
    """
    Recibe un ángulo theta y retorna una matriz de 2x2
    que rota un vector dado en un ángulo theta
    """
    res = np.zeros((2,2))
    res[0][0]=np.cos(theta)
    res[0][1]=-(np.sin(theta))
    res[1][0]=np.sin(theta)
    res[1][1]=np.cos(theta)
    return res


def escala(s):
    """
    Recibe una tira de números s y retorna una matriz cuadrada de
    n x n, donde n es el tamaño de s.
    La matriz escala la componente i de un vector de Rn
    en un factor s[i]
    """
    res = np.zeros((len(s),len(s)))
    for i in range(len(s)):
        res[i][i]=s[i]
    return res


def rota_y_escala(theta, s):
    """
    Recibe un ángulo theta y una tira de números s,
    y retorna una matriz de 2x2 que rota el vector en un ángulo theta
    y luego lo escala en un factor s
    """
    res=multiplicar_matrices(escala(s),rota(theta))
    return np.array(res)


def afin(theta, s, b):
    """
    Recibe un ángulo theta, una tira de números s (en R2), y un vector b en R2.
    Retorna una matriz de 3x3 que rota el vector en un ángulo theta,
    luego lo escala en un factor s y por último lo mueve en un valor fijo b
    """
    rota_escala=rota_y_escala(theta,s)
    res=np.zeros((3,3))
    res[0][0]=rota_escala[0][0]
    res[0][1]=rota_escala[0][1]
    res[1][0]=rota_escala[1][0]
    res[1][1]=rota_escala[1][1]
    res[0][2]=b[0]
    res[1][2]=b[1]
    res[2][2]=1
    return res


def trans_afin(v, theta, s, b):
    """
    Recibe un vector v (en R2), un ángulo theta,
    una tira de números s (en R2), y un vector b en R2.
    Retorna el vector w resultante de aplicar la transformación afín a v
    """
    SRT=afin(theta,s,b)
    v_columna=[v[0],v[1],v[2]]
    w_columna=multiplicar_matrices(SRT,v_columna)
    return [w_columna[0][0], w_columna[1][0]]

#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

# Labo 3

def norma(x, p):
    """
    Devuelve la norma p del vector x.
    """
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
    """
    Recibe X, una lista de vectores no vacíos, y un escalar p. Devuelve
    una lista donde cada elemento corresponde a normalizar los
    elementos de X con la norma p.
    """
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


def normaMatMC(A, q, p, Np):
    """
    Devuelve la norma ||A||_{q,p} y el vector x en el cual se alcanza
    el máximo.
    """
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



def normaExacta(A, p=[1, 'inf']):
    """
    Devuelve una lista con las normas 1 e infinito de una matriz A
    usando las expresiones del enunciado 2.(c).
    """
        """Calcula la norma matricial subordinada para p=1 o p='inf'."""
    if A.ndim != 2 or A.size == 0:
        return None

    if p == 1:
        return norma1(A)

    if p == 'inf':
        return normaInf(A)

    return None



def condMC(A, p):
    """
    Devuelve el número de condición de A usando la norma inducida p.
    """
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
    """
    Que devuelve el número de condición de A a partir de la fórmula de
    la ecuación (1) usando la norma p.
    """
    A_inv = inversa(A)

    if A_inv is None:
        return None

    norma_A = normaExacta(A, p)
    norma_Ainv = normaExacta(A_inv, p)

    if norma_A is None or norma_Ainv is None:
        return None

    return float(norma_A * norma_Ainv)


# Labo 4

def calculaLU(A):
    """
    Calcula la factorización LU de la matriz A y retorna las matrices L
    y U, junto con el número de operaciones realizadas. En caso de
    que la matriz no pueda factorizarse retorna None.
    """
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

def restri(L, b, inferior=True):
    """
    Resuelve el sistema Lx = b, donde L es triangular. Se puede indicar
    si es triangular inferior o superior usando el argumento
    inferior (por default asumir que es triangular inferior).
    """
   n = b.shape[0]

    x = np.zeros(n, dtype=float)

    if inferior:
        for i in range(n):
            if abs(T[i, i]) <= 1e-15:
                return None

            suma = 0.0
            for j in range(i):
                suma += T[i, j] * x[j]

            x[i] = (b[i] - suma) / T[i, i]
    else:
        for i in range(n - 1, -1, -1):
            if abs(T[i, i]) <= 1e-15:
                return None

            suma = 0.0
            for j in range(i + 1, n):
                suma += T[i, j] * x[j]

            x[i] = (b[i] - suma) / T[i, i]

    return x


def inversa(A):
    """
    Calcula la inversa de A empleando la factorización LU
    y las funciones que resuelven sistemas triangulares.
    """
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
    """
    Calcula la factorización LDV de la matriz A, de forma tal que A =
    LDV, con L triangular inferior, D diagonal y V triangular
    superior. En caso de que la matriz no pueda factorizarse
    retorna None.
    """
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
    """
    Checkea si la matriz A es simétrica definida positiva (SDP) usando
    la factorización LDV.
    """
    if A.ndim != 2:
        return False

    n, m = A.shape
    if n != m or n == 0:
        return False

    # Verificacion explicita de simetria
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i, j] - A[j, i]) > atol:
                return False

    # Obtencion de factorizacion LDV
    res = calculaLDV(A)
    if res is None:
        return False

    L, D, V = res

    # Verificacion de elementos diagonales positivos
    for i in range(n):
        if D[i, i] <= atol:
            return False

    return True


#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#

# Labo 5

def QR_con_GS(A, tol=1e-12, retorna_nops=False):
    """
    A una matriz de n x n
    tol la tolerancia con la que se filtran elementos nulos en R
    retorna_nops permite (opcionalmente) retornar el numero de
    operaciones realizado
    retorna matrices Q y R calculadas con Gram Schmidt (y como
    tercer argumento opcional, el numero de operaciones).
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


def QR_con_HH(A, tol=1e-12, extras=False):
    """
    A una matriz de m x n (m>=n)
    tol la tolerancia con la que se filtran elementos nulos en R
    retorna matrices Q y R calculadas con reflexiones de
    Householder
    Si la matriz A no cumple m>=n, debe retornar None
    extras : bool, opcional
    Si es True, devuelve informacion extra sobre el proceso de
    factorizacion.
    Por defecto es False. Esto lo hacemos para poder graficar
    el proceso.
    Devuelve la factorizacion QR de A usando reflectores de
    Householder.
    Devuelve:
    Q, R, extra_info (si extras es True)
    Q, R (si extras es False)
    extra_info es un diccionario con la clave:
    'R_matrices' : lista de las matrices R en cada paso
    'Q_matrices' : lista de las matrices Q en cada paso
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



def calculaQR(A, metodo='RH', tol=1e-12):
    """
    A una matriz de n x n
    tol la tolerancia con la que se filtran elementos nulos en R
    metodo = ['RH', 'GS'] usa reflectores de Householder (RH) o Gram
    Schmidt (GS) para realizar la factorizacion
    retorna matrices Q y R calculadas con Gram Schmidt (y como
    tercer argumento opcional, el numero de operaciones)
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


# Labo 6

def metpot2k(A, tol=1e-15, K=1000):
    """
    A una matriz de n x n
    tol la tolerancia en la diferencia entre un paso y el siguiente de la estimacion del
    autovector.
    K el numero maximo de iteraciones a realizarse.
    Retorna vector v, autovalor lambda y numero de iteracion realizadas k
    """
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

def diagRH(A, tol=1e-15, K=1000):
    """
    A una matriz simetrica de n x n
    tol la tolerancia en la diferencia entre un paso y el siguiente de la estimacion del
    autovector. K el numero maximo de iteraciones a realizarse
    retorna matriz de autovectores S y matriz de autovalores D, tal que A = S D S.T
    Si la matriz A no es simetrica, debe retornar None
    """
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

# Labo 7

def transicionesalazarcontinuas(n):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    Retorna matriz T de n x n normalizada por columnas, y con entradas al azar en el
    intervalo [0, 1]
    """
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            res[i][j] = random.uniform(0.0, 1.0)

    return res / np.sum(res, axis=0)   

def transicionesalazaruniforme(n, thres):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    thres probabilidad de que una entrada sea distinta de cero.
    Retorna matriz T de n x n normalizada por columnas. El elemento i, j es distinto de
    cero si el número generado al azar para i, j es menor o igual a thres. Todos los
    elementos de la columna $j$ son iguales (a 1 sobre el número de elementos distintos
    de cero en la columna).
    """
    res = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            res[i][j] = random.uniform(0.0, 1.0)

    for i in range(n):
        for j in range(n):
            if res[i][j] <= thres:
                res[i][j] = 1
            else: res[i][j] = 0

    for j in range(n):
        suma_col = 0
        for i in range(n):
            suma_col += res[i][j]
        if suma_col == 0:
            indice = random.randint(0, n-1)
            res[indice][j] = 1.0    

    return res / np.sum(res, axis=0)


def nucleo(A, tol=1e-15):
    """
    A una matriz de m x n
    tol la tolerancia para asumir que un vector esta en el nucleo.
    Calcula el nucleo de la matriz A diagonalizando la matriz traspuesta(A) * A (* la
    multiplicacion matricial), usando el medodo diagRH. El nucleo corresponde a los
    autovectores de autovalor con modulo <= tol.
    Retorna los autovectores en cuestion, como una matriz de n x k, con k el numero de
    autovectores en el nucleo.
    """
    m, n = A.shape
    B = arroba(A.T, A) # B de dim nxn
    S,D = diagRH(B, tol)

    indices = []
    k = 0
    for i in range(n):
        if np.abs(D[i][i]) <= tol:
            indices.append(i)
            k += 1

    if k == 0:
        return np.zeros((0, 0))

    res = np.zeros((n, k))
    for j in range(k):
        l = indices[j]
        for i in range(n):
            res[i][j] = S[i][l]
    return res


def crearabla(listado, mfilas, ncolumnas, tol=1e-15):
    """
    Recibe una lista listado, con tres elementos: lista con indices i, lista con indices
    j, y lista con valores Aij de la matriz A. Tambien las dimensiones de la matriz a
    traves de mfilas y ncolumnas. Los elementos menores a tol se descartan.
    Idealmente, el listado debe incluir unicamente posiciones correspondientes a valores
    distintos de cero. Retorna una lista con:
    - Diccionario {(i, j) : Aij} que representa los elementos no nulos de la matriz A. Los
    elementos con modulo menor a tol deben descartarse por default.
    - Tupla (mfilas, ncolumnas) que permita conocer las dimensiones de la matriz.
    """
    dim = (m_filas, n_columnas)
    dic = dict()

    if len(listado) == 0 or len(listado[0]) == 0: 
        return [dic, dim] 
    
    indices_i = listado[0]
    indices_j = listado[1]
    valores_Aij = listado[2]
    for k in range(len(indices_i)):
        if abs(valores_Aij[k]) > tol:
            dic[indices_i[k], indices_j[k]] = valores_Aij[k]
    return [dic, dim]

def multiplicaralavector(A, v):
    """
    Recibe una matriz rala creada con crearabla y un vector v.
    Retorna un vector w resultado de multiplicar A con v
    """
    m, n = A[1]
    v = np.asarray(v)

    if n != len(v):
        return

    res = np.zeros(m)
    for i,j in A[0]:
        res[i] += A[0][i,j] * v[j]

    return res

# Labo 8

def svdreducida(A, k="max", tol=1e-15):
    """
    A la matriz de interes (de m x n)
    k el numero de valores singulares (y vectores) a retener.
    tol la tolerancia para considerar un valor singular igual a cero
    Retorna hatU (matriz de m x k), hatSig (vector de k valores singulares) y hatV (
    matriz de n x k)
    """
    m, n = A.shape
    
    # Obtenemos los componentes base de la función auxiliar
    res1, valores_singulares, res3, r, transpuesta = _calcular_bases_svd(A, tol)
    
    # Definir el número efectivo de componentes k a retener
    if k == "max":
        k = max(1, r)
    else:
        k = min(int(k), len(valores_singulares), r)
        
    # Recortar los valores singulares a k dimensiones
    hatSig = valores_singulares[:k]
    
    # Inicializar las matrices finales truncadas
    hatU = np.zeros((m, k))
    hatV = np.zeros((n, k))
    
    if not transpuesta:
        # Caso m >= n: 
        # res1 es V (n x n), res3 es B = U * Sigma (m x n)
        V = res1
        U_sig = res3
        
        hatV = V[:, :k]
        for i in range(k):
            if hatSig[i] > tol:
                hatU[:, i] = U_sig[:, i] / hatSig[i]
            else:
                hatU[:, i] = 0.0
    else:
        # Caso m < n (estaba transpuesta): 
        # res1 es B = V * Sigma (n x m), res3 es U (m x m)
        V_sig = res1
        U = res3
        
        hatU = U[:, :k]
        for i in range(k):
            if hatSig[i] > tol:
                hatV[:, i] = V_sig[:, i] / hatSig[i]
            else:
                hatV[:, i] = 0.0
                
    return hatU, hatSig, hatV



#-----------------FUNCIONES AUXILIARES---------------

def multiplicar_matrices(A, B):
    filas_A = len(A)
    columnas_A = len(A[0])
    filas_B = len(B)
    columnas_B = len(B[0])

    resultado = [[0 for _ in range(columnas_B)] for _ in range(filas_A)]

    for i in range(filas_A):
        for j in range(columnas_B):
            for k in range(columnas_A):
                resultado[i][j] += A[i][k] * B[k][j]

    return resultado


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


def householder(v1):
    n = len(v1)
    e1 = np.zeros(n)
    e1[0] = 1.0
    u = e1 - v1
    nu2 = u @ u
    if nu2 < 1e-30:                  # v1 ya es e1: H = I (evita dividir por 0)
        return np.eye(n)
    return np.eye(n) - 2.0 * np.outer(u, u) / nu2


def arroba(A, B):
    # Multiplicación matricial sin @
    A = np.asarray(A)
    B = np.asarray(B)
    
    if A.ndim == 2 and B.ndim == 2:
        res = np.zeros((A.shape[0], B.shape[1]))
        for i in range(A.shape[0]):
            for j in range(B.shape[1]):
                res[i, j] = np.sum(A[i, :] * B[:, j])
        return res
    elif A.ndim == 2 and B.ndim == 1:
        res = np.zeros(A.shape[0])
        for i in range(A.shape[0]):
            res[i] = np.sum(A[i, :] * B)
        return res
    elif A.ndim == 1 and B.ndim == 1:
        return np.sum(A * B)


def _calcular_bases_svd(A, tol=1e-15):
    """
    Función interna que evalúa dimensiones, diagonaliza y calcula 
    los componentes base de la SVD (valores singulares, V y B).
    """
    m, n = A.shape
    transpuesta = False
    
    # Optimización por dimensiones (m < n vs m > n)
    if m < n:
        A = A.T
        transpuesta = True
        
    # Diagonalización de A^t * A
    AtA = np.dot(A.T, A) #producto entre A y A^t
    autovalores, autovectores = diagRH(AtA)  # obtengo autovalores y autovectores de A * At
    
    # Ordenar de mayor a menor
    idx_ordenados = np.argsort(autovalores)[::-1] 
    autovalores = autovalores[idx_ordenados]
    autovectores = autovectores[:, idx_ordenados]
    
    # Calcular valores singulares
    valores_singulares = np.sqrt(np.maximum(autovalores, 0))
    
    # Rango efectivo (valores mayores a la tolerancia)
    r = np.sum(valores_singulares > tol)
    
    # Calcular B = A * V para luego obtener U
    B = np.dot(A, autovectores)
    
    if transpuesta:
        # Si se traspuso, intercambiamos roles para devolver U y V correctamente ordenados
        return B, valores_singulares, autovectores, r, True
    else:
        return autovectores, valores_singulares, B, r, False
