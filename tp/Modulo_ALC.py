# Labo 0 



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

# Labo 5


# Labo 6


# Labo 7


# Labo 8


