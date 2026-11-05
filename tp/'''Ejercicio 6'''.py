'''Ejercicio 6'''


import numpy as np
import alc as alc


def productoCruz(a, b):
    """Producto vectorial de dos vectores de R^3."""
    c = np.zeros(3)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c


def calculoLambdaMu(I0, J0):
    """
    Resuelve lambda y mu de la ecuación 16 (C^2 = (lambda + i mu)^2).
    Parámetros: I0, J0 (vectores de R^3), soluciones del paso POS.
    Retorna: (lam, mu) con lam >= 0. La otra solución es (-lam, -mu).
    Se cumple: lam*mu = -I0.J0   y   lam^2 - mu^2 = |J0|^2 - |I0|^2.
    """

    '''Dos restricciones I^T J = 0 y ||I||^2 = ||J||^2.
    lambda ^2 - u^2 = ||J_0||^2 - ||I_0||^2
    C = lambda = i.u ====> C^2 = lamda ^2 - U^2 + 2i lamda U ====> C^2 = ||J_0||^2 - ||I_0||^2 -2i I_0 J_0'''
    I0 = np.array(I0)
    J0 = np.array(J0)
    r = np.sum(J0 * J0) - np.sum(I0 * I0)       # parte real de C^2 = ||J_0||^2 - ||I_0||^2
    i = -2.0 * np.sum(I0 * J0)                  # parte imaginaria de C^2 = -2i I_0 J_0
    rho = np.sqrt(r * r + i * i)             # módulo de C^2
    lam = np.sqrt((rho + r) / 2.0)
    if lam > 1e-12:
        mu = i / (2.0 * lam)
    else:                                         # caso lam = 0: C^2 real y <= 0
        mu = np.sqrt(max((rho - r) / 2.0, 0.0))
    return lam, mu


def calcularSalidaPOSIT(Is, Js, fx, m0, M0):
    """
    Pose asociada a una solución (Is, Js) del paso POS.
    Parámetros:
        Is, Js: vectores I y J (R^3).
        fx: distancia focal en píxeles.
        m0: (x0, y0) imagen de M0 centrada en el punto principal.
        M0: coordenadas de M0 en el sistema del objeto (R^3).
    Retorna: i, j, k, Z0, R, t   (X_cam = R*M + t).
    """
    Is = np.array(Is)
    Js = np.array(Js)
    nI = alc.norma(Is, 2)
    nJ = alc.norma(Js, 2)
    Z0 = fx / nI                                  # = fx / ||J|| si las normas coinciden
    i = Is / nI
    j = Js / nJ
    k = productoCruz(i, j)
    R = np.zeros((3, 3))
    R[0, :] = i
    R[1, :] = j
    R[2, :] = k
    # T = posición de M0 en coordenadas de cámara
    T = np.array([Z0 * m0[0] / fx, Z0 * m0[1] / fx, Z0])
    # X_cam = R (M - M0) + T  =>  t = T - R M0
    t = T - alc.arroba(R, np.asarray(M0))
    return i, j, k, Z0, R, t


def evaluarPose(R, t, obj_vecs, img_vecs, K):
    """
    Cuenta las proyecciones válidas (Z > 0) y calcula el RMS de reproyección
    (en píxeles) sobre los puntos válidos.
    Retorna: (cantidad_validas, rms). Si no hay válidas, rms = inf.
    """
    obj_vecs = np.array(obj_vecs)
    img_vecs = np.array(img_vecs)
    n = obj_vecs.shape[0]
    validas = 0
    suma = 0.0
    for p in range(n):
        Xc = alc.arroba(R, obj_vecs[p]) + t
        if Xc[2] > 0:
            validas = validas + 1
            u = (K[0, 0] * Xc[0] + K[0, 1] * Xc[1]) / Xc[2] + K[0, 2]
            v = K[1, 1] * Xc[1] / Xc[2] + K[1, 2]
            du = u - img_vecs[p, 0]
            dv = v - img_vecs[p, 1]
            suma = suma + du * du + dv * dv
    if validas == 0:
        return 0, np.inf
    return validas, np.sqrt(suma / validas)


def pseudoinversa_rango2(A):
    """Pseudoinversa A^+ = V S^-1 U^T de una matriz de rango 2, y la normal al plano."""
    U, S, V = alc.svd_reducida(A, k=2)
    n_ = V.shape[0]
    Ap = np.zeros((n_, A.shape[0]))
    for c in range(2):
        Ap = Ap + (V[:, c].reshape(n_, 1) * U[:, c].reshape(1, A.shape[0])) / S[c]
    normal = productoCruz(V[:, 0], V[:, 1])
    normal = normal / alc.norma(normal, 2)
    return Ap, normal


def aplicarPOSIT(obj_vecs, img_vecs, K, MaxIt=50, atol=1e-6):
    """
    POSIT coplanar. obj_vecs (n x 3): puntos del objeto; el primero es M0.
    img_vecs (n x 2): píxeles (u, v) de los mismos puntos, mismo orden.
    Retorna (R, t) con X_cam = R*M + t (t en las unidades de obj_vecs), o None.
    """
    obj_vecs = np.asarray(obj_vecs, dtype=float)
    img_vecs = np.asarray(img_vecs, dtype=float)
    n = obj_vecs.shape[0]
    fx = K[0, 0]

    # Coordenadas de imagen centradas; y se reescala por fx/fy para tener una
    # cámara equivalente con una única focal fx
    x = img_vecs[:, 0] - K[0, 2]
    y = (img_vecs[:, 1] - K[1, 2]) * (fx / K[1, 1])
    x0, y0 = x[0], y[0]
    M0 = obj_vecs[0]

    A = obj_vecs[1:] - M0                          # a_i = M0Mi
    Ap, normal = pseudoinversa_rango2(A)

    eps = np.zeros(n - 1)
    rms_ant = np.inf
    mejor = None
    for it in range(MaxIt):
        bu = x[1:] * (1.0 + eps) - x0
        bv = y[1:] * (1.0 + eps) - y0
        I0 = alc.arroba(Ap, bu)
        J0 = alc.arroba(Ap, bv)
        lam, mu = calculoLambdaMu(I0, J0)

        elegida = None
        for signo in (1.0, -1.0):                  # las dos poses simétricas
            Is = I0 + signo * lam * normal
            Js = J0 + signo * mu * normal
            i, j, k, Z0, R, t = calcularSalidaPOSIT(Is, Js, fx, (x0, y0), M0)
            validas, rms = evaluarPose(R, t, obj_vecs, img_vecs, K)
            if validas == n and (elegida is None or rms < elegida[0]):
                elegida = (rms, R, t, k, Z0)
        if elegida is None:
            break                                  # ninguna pose con todos los puntos delante
        rms, R, t, k, Z0 = elegida
        mejor = (R, t)
        eps = alc.arroba(A, k) / Z0                # nuevos factores de corrección
        if abs(rms_ant - rms) < atol:
            break
        rms_ant = rms
    return mejor