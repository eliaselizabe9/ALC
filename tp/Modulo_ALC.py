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

def error relativo(x,y):
    '''
    Recibe dos numeros x e y, y calcula el error relativo de aproximar x usando y en float64
    '''

def matricesIguales(A,B):
    '''
    Devuelve True si ambas matrices son iguales y False en otro caso.
    Considerar que las matrices pueden tener distintas dimensiones, ademas de distintos valores.
    '''

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

def escala(s):
    """
    Recibe una tira de números s y retorna una matriz cuadrada de
    n x n, donde n es el tamaño de s.
    La matriz escala la componente i de un vector de Rn
    en un factor s[i]
    """

def rota_y_escala(theta, s):
    """
    Recibe un ángulo theta y una tira de números s,
    y retorna una matriz de 2x2 que rota el vector en un ángulo theta
    y luego lo escala en un factor s
    """

def afin(theta, s, b):
    """
    Recibe un ángulo theta, una tira de números s (en R2), y un vector b en R2.
    Retorna una matriz de 3x3 que rota el vector en un ángulo theta,
    luego lo escala en un factor s y por último lo mueve en un valor fijo b
    """

def transafin(v, theta, s, b):
    """
    Recibe un vector v (en R2), un ángulo theta,
    una tira de números s (en R2), y un vector b en R2.
    Retorna el vector w resultante de aplicar la transformación afín a v
    """

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

def normaliza(X, p):
    """
    Recibe X, una lista de vectores no vacíos, y un escalar p. Devuelve
    una lista donde cada elemento corresponde a normalizar los
    elementos de X con la norma p.
    """

def normaMatMC(A, q, p, Np):
    """
    Devuelve la norma ||A||_{q,p} y el vector x en el cual se alcanza
    el máximo.
    """

def normaExacta(A, p=[1, 'inf']):
    """
    Devuelve una lista con las normas 1 e infinito de una matriz A
    usando las expresiones del enunciado 2.(c).
    """

def condMC(A, p):
    """
    Devuelve el número de condición de A usando la norma inducida p.
    """

def condExacto(A, p):
    """
    Que devuelve el número de condición de A a partir de la fórmula de
    la ecuación (1) usando la norma p.
    """

# Labo 4

def calculaLU(A):
    """
    Calcula la factorización LU de la matriz A y retorna las matrices L
    y U, junto con el número de operaciones realizadas. En caso de
    que la matriz no pueda factorizarse retorna None.
    """

def restri(L, b, inferior=True):
    """
    Resuelve el sistema Lx = b, donde L es triangular. Se puede indicar
    si es triangular inferior o superior usando el argumento
    inferior (por default asumir que es triangular inferior).
    """

def inversa(A):
    """
    Calcula la inversa de A empleando la factorización LU
    y las funciones que resuelven sistemas triangulares.
    """

def calculaLDV(A):
    """
    Calcula la factorización LDV de la matriz A, de forma tal que A =
    LDV, con L triangular inferior, D diagonal y V triangular
    superior. En caso de que la matriz no pueda factorizarse
    retorna None.
    """

def esSDP(A, atol=1e-8):
    """
    Checkea si la matriz A es simétrica definida positiva (SDP) usando
    la factorización LDV.
    """

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

# Labo 6

def metpot2k(A, tol=1e-15, K=1000):
    """
    A una matriz de n x n
    tol la tolerancia en la diferencia entre un paso y el siguiente de la estimacion del
    autovector.
    K el numero maximo de iteraciones a realizarse.
    Retorna vector v, autovalor lambda y numero de iteracion realizadas k
    """

def diagRH(A, tol=1e-15, K=1000):
    """
    A una matriz simetrica de n x n
    tol la tolerancia en la diferencia entre un paso y el siguiente de la estimacion del
    autovector. K el numero maximo de iteraciones a realizarse
    retorna matriz de autovectores S y matriz de autovalores D, tal que A = S D S.T
    Si la matriz A no es simetrica, debe retornar None
    """


# Labo 7

def transicionesalazarcontinuas(n):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    Retorna matriz T de n x n normalizada por columnas, y con entradas al azar en el
    intervalo [0, 1]
    """

def transicionesalazaruniforme(n, thres):
    """
    n la cantidad de filas (columnas) de la matriz de transición.
    thres probabilidad de que una entrada sea distinta de cero.
    Retorna matriz T de n x n normalizada por columnas. El elemento i, j es distinto de
    cero si el número generado al azar para i, j es menor o igual a thres. Todos los
    elementos de la columna $j$ son iguales (a 1 sobre el número de elementos distintos
    de cero en la columna).
    """

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

def multiplicaralavector(A, v):
    """
    Recibe una matriz rala creada con crearabla y un vector v.
    Retorna un vector w resultado de multiplicar A con v
    """


# Labo 8

def svdreducida(A, k="max", tol=1e-15):
    """
    A la matriz de interes (de m x n)
    k el numero de valores singulares (y vectores) a retener.
    tol la tolerancia para considerar un valor singular igual a cero
    Retorna hatU (matriz de m x k), hatSig (vector de k valores singulares) y hatV (
    matriz de n x k)
    """

