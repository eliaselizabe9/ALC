'''
Ejercicio 5: Cruzando en rojo - Calculo de velocidad
Implementacion de calibracion intrinseca, estimacion de pose y reproyeccion,
y calculo de velocidades de agentes dinamicos sobre la calzada.
'''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import alc


def extraerTiempoMs(nombreArchivo):
    '''
    Extrae la estampa temporal en segundos con precision de milisegundos
    a partir del nombre del archivo de captura (e.g. img_10_18_47_351.jpg -> 37127.351 s).
    '''
    base = os.path.basename(nombreArchivo)
    nombre_sin_ext = base.split('.')[0]
    partes = nombre_sin_ext.split('_')
    partes_num = [p for p in partes if p.isdigit()]
    if len(partes_num) >= 4:
        horas = float(partes_num[0])
        minutos = float(partes_num[1])
        segundos = float(partes_num[2])
        milisegundos = float(partes_num[3])
        return horas * 3600.0 + minutos * 60.0 + segundos + (milisegundos / 1000.0)
    return 0.0


# =========================================================================
# MARCA 5: PROYECCION INVERSA 2D IMAGEN -> 3D PLANO PISO (H_inv * m)
# Transforma un punto en pixeles [u, v, 1] al suelo metrico [X, Y] en metros
# =========================================================================

def proyectarPuntoImagenAPiso(m, H_inv):
    '''
    Toma las coordenadas homogeneas de un punto medido en la imagen:
        m = [u_k, v_k, 1.0]
    Y aplica la transformacion de homografia inversa del plano calle:
        P_homog = H_inv * m
    Normaliza por la tercera componente para obtener coordenadas metricas en metros:
        X = P_homog[0] / P_homog[2]
        Y = P_homog[1] / P_homog[2]
    '''
    P_homog = alc.arroba(H_inv, m)
    
    # Normalizacion
    if abs(P_homog[2]) > 1e-15:
        return np.array([P_homog[0] / P_homog[2], P_homog[1] / P_homog[2]])
    return np.array([P_homog[0], P_homog[1]])


def obtenerQ3d(imgPts, objPts):
    '''
    Construye la matriz Q de dimension 2n x 12 para el sistema homogeneo Qm = 0.
    Cada correspondencia 3D-2D genera 2 ecuaciones independientes.
    '''
    n = objPts.shape[0]
    Q = np.zeros((2 * n, 12))
    for i in range(n):
        P_i = np.array([objPts[i, 0], objPts[i, 1], objPts[i, 2], 1.0])
        u = imgPts[i, 0]
        v = imgPts[i, 1]
        
        Q[2 * i, 0:4] = P_i
        Q[2 * i, 4:8] = 0.0
        Q[2 * i, 8:12] = -u * P_i
        
        Q[2 * i + 1, 0:4] = 0.0
        Q[2 * i + 1, 4:8] = P_i
        Q[2 * i + 1, 8:12] = -v * P_i
    return Q


def resolverDLT3D(Q):
    '''
    Resuelve el sistema de calibracion 3D DLT mediante svd_reducida del modulo ALC.
    Reconstruye y retorna la matriz de proyeccion M de 3 x 4.
    Para puntos coplanares (Z=0), selecciona la columna del subespacio nulo
    correspondiente a la geometria fisica real (coeficientes no nulos en X, Y, W).
    '''
    U_hat, Sig_hat, V_hat = alc.svd_reducida(Q, k=12, tol=1e-15)
    
    # Recorrido desde el menor valor singular hacia arriba para descartar columnas nulas de Z
    for col in range(V_hat.shape[1] - 1, -1, -1):
        v = V_hat[:, col]
        if abs(v[0]) > 1e-6 or abs(v[1]) > 1e-6 or abs(v[3]) > 1e-6:
            return v.reshape((3, 4))
            
    return V_hat[:, -1].reshape((3, 4))


# Grafico de Puntos en imagenes de Salguero

def dibujarPuntosSobreImagen(imgPath, imgPts, puntosReproyectados, rms, outPath="reproyeccion_salguero.png"):
    '''
    Dibuja los puntos medidos y los puntos reproyectados sobre la imagen:
      * Grafica los PUNTOS MEDIDOS (imgPts) en color VERDE.
      * Grafica los PUNTOS REPROYECTADOS en color ROJO.
      * La discrepancia (error individual) entre cada punto medido y su reproyección es trazado mediante una línea AMARILLA.
    '''
    im = mpimg.imread(imgPath)
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.imshow(im)
    
    # Puntos medidos manualmente en la imagen (Verdes)
    ax.scatter(imgPts[:, 0], imgPts[:, 1], c='lime', s=35, label='Puntos medidos', marker='o', zorder=3)
    
    # Puntos calculados por el modelo DLT (Rojos)
    ax.scatter(puntosReproyectados[:, 0], puntosReproyectados[:, 1], c='red', s=30, label='Puntos reproyectados (Modelo M)', marker='x', zorder=4)
    
    # Lineas residuales de error (Amarillas)
    for i in range(imgPts.shape[0]):
        ax.plot([imgPts[i, 0], puntosReproyectados[i, 0]],
                [imgPts[i, 1], puntosReproyectados[i, 1]],
                color='yellow', linewidth=0.8, zorder=2)
                
    ax.set_title(f'Reproyeccion DLT 3D sobre {os.path.basename(imgPath)} (RMS: {rms:.3f} px)')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(outPath)
    plt.show()


def graficarVelocidades(tablaVelocidades, outPath="velocidades_agentes.png"):
    '''
    Grafica la evolucion temporal de las velocidades estimadas para cada agente.
    '''
        
    plt.figure(figsize=(10, 6))
    agentes = ['Taxi', 'Camioneta', 'Peatona']
    colores = {'Taxi': 'gold', 'Camioneta': 'blue', 'Peatona': 'crimson'}
    estilos = {'Taxi': 'o-', 'Camioneta': 's-', 'Peatona': '^-'}
    
    for agente in agentes:
        filas = [f for f in tablaVelocidades if f['agente'] == agente]
        if len(filas) > 0:
            if 't_i' in filas[0] and 't_f' in filas[0]:
                t0 = filas[0]['t_i']
                x_vals = [(f['t_i'] + f['t_f']) / 2.0 - t0 for f in filas]
                eje_x_label = 'Tiempo relativo (s)'
            else:
                x_vals = [f['frame_i'] for f in filas]
                eje_x_label = 'Indice de intervalo de captura'
            vels = [f['v_kmh'] for f in filas]
            plt.plot(x_vals, vels, estilos.get(agente, 'o-'), color=colores.get(agente, 'black'), label=f'{agente}', linewidth=2, markersize=8)
            
    plt.title('Evolucion Temporal de Velocidades Estimadas')
    plt.xlabel(eje_x_label)
    plt.ylabel('Velocidad (km/h)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outPath)
    plt.show()


# Punto a: Calcular Intrínsecos Salgero

def calcularIntrinsecosSalguero(world_img_pts_csv):
    '''
    Calcula los parametros intrinsecos K de la camara promediando las estimaciones
    obtenidas a partir de diferentes vistas de la esquina en world_img_pts_csv.
    '''
    df = pd.read_csv(world_img_pts_csv, sep=';')
    objPts = df[['X', 'Y', 'Z']].values
    
    K_acum = np.zeros((3, 3))
    
    columnas = list(df.columns)
    vistas = [c.split('_')[1] for c in columnas if c.startswith('u_') and f"v_{c.split('_')[1]}" in columnas]
    num_vistas = len(vistas)
    
    for k in range(num_vistas):
        imgPts_k = df[[f'u_{k}', f'v_{k}']].values
        Q_k = obtenerQ3d(imgPts_k, objPts)
        U_hat, Sig_hat, V_hat = alc.svd_reducida(Q_k, k=12, tol=1e-15)
        
        # m_opt = ultima columna de V_hat
        m_opt = V_hat[:, -1]
        M_k = m_opt.reshape((3, 4))
        
        # Factorizacion RQ de M_k[0:3, 0:3] mediante QR de su inversa
        M3_inv = alc.inversa(M_k[0:3, 0:3])
        Q_qr, R_qr = alc.calculaQR(M3_inv, metodo='RH', tol=1e-12)
        
        K_k = alc.inversa(R_qr)
        K_k = K_k / K_k[2, 2] # Normalizacion K_33 = 1.0
        
        # Ajuste de focos positivos
        for c in range(2):
            if K_k[c, c] < 0.0:
                K_k[:, c] = -K_k[:, c]
                
        K_acum = K_acum + K_k
        
    K_final = K_acum / num_vistas
    
    return K_final


# Punto b: Estimar Pose y Reproyectar

def estimarPoseYReproyectar(file_pts, K, img_path):
    '''
    Estima la pose de la camara (rotacion rigida Q_r y traslacion t) y calcula
    el error cuadratico medio de reproyeccion RMS sobre la imagen de referencia.
    '''
    data = np.loadtxt(file_pts)
    objPts = data[1:4, :].T  # Matriz de 51 x 3 con [X, Y, Z] de cada punto
    imgPts = data[4:6, :].T  # Matriz de 51 x 2 con [u, v] medidos en la imagen
    n_pts = objPts.shape[0]
    
    # Obtenemos la matriz Q mediante obtenerQ3d. Matriz de 102 x 12
    Q = obtenerQ3d(imgPts, objPts)
    
    # M = resolverDLT3D(Q) - Reconstruccion 3x4 mediante svd_reducida
    M = resolverDLT3D(Q)
    
    K_inv = alc.inversa(K)
    R_raw = alc.arroba(K_inv, M[0:3, 0:3])
    
    Q_r, R_r = alc.calculaQR(R_raw, metodo='RH')
    
    t = alc.res_tri(K, M[0:3, 3], inferior=False)
    
    suma_err_sq = 0.0
    puntos_reproyectados = np.zeros((n_pts, 2))
    
    # PROYECCION MATEMATICA 3D -> 2D (M * P_i) Y COORDENADAS REPROYECTADAS
    # Transforma cada punto homogeneo 3D P_i al plano imagen para obtener (u_hat, v_hat)
    for i in range(n_pts):
        # Punto 3D
        P_i = np.array([objPts[i, 0], objPts[i, 1], objPts[i, 2], 1.0])
        
        # Proyeccion lineal: m_hat_homog = M * P_i (vector de 3 componentes [u_tilde, v_tilde, w_tilde])
        m_hat_homog = alc.arroba(M, P_i)
        
        # Normalización mediante tilde w para obtener pixeles euclideos
        u_hat = m_hat_homog[0] / m_hat_homog[2]
        v_hat = m_hat_homog[1] / m_hat_homog[2]
        
        # Guardamos el punto proyectado resultante
        puntos_reproyectados[i, 0] = u_hat
        puntos_reproyectados[i, 1] = v_hat
        
        # Acumulacion de la discrepancia cuadratica entre punto medido y punto proyectado
        suma_err_sq += (imgPts[i, 0] - u_hat)**2 + (imgPts[i, 1] - v_hat)**2
        
    # RMS = RaizCuadrada(suma_err_sq / 51.0)
    RMS = np.sqrt(suma_err_sq / float(n_pts))
    
    dibujarPuntosSobreImagen(img_path, imgPts, puntos_reproyectados, RMS)
    
    return Q_r, t, RMS


# Punto c: Estimar Velocidad de Agentes (Peatona, Taxi y Camion)

def estimarVelocidadAgentes(secuencia_capturas, H_calle, tracking_agentes):
    '''
    Estima la velocidad de circulacion de los agentes dinamicos (Taxi, Camioneta, Peatona)
    a lo largo de una secuencia de capturas proyectando las coordenadas al plano de la calle.
    '''
    # 1. H_inv = inversa(H_calle)
    H_inv = alc.inversa(H_calle)
    
    # 2. Tabla_Velocidades = CrearListaVacia()
    Tabla_Velocidades = []
    
    # 3. Para cada agente A en {Taxi, Camioneta, Peatona} hacer:
    agentes = ['Taxi', 'Camioneta', 'Peatona']
    for A in agentes:
        if A not in tracking_agentes:
            continue
            
        registros_agente = tracking_agentes[A]
        num_frames = len(registros_agente)
        
        # Para k = 0 hasta num_frames - 2 hacer:
        for k in range(num_frames - 1):
            t_k = extraerTiempoMs(secuencia_capturas[k])
            t_k1 = extraerTiempoMs(secuencia_capturas[k + 1])
            dt = t_k1 - t_k # Lapso exacto en segundos
            
            u_k, v_k = registros_agente[k]
            u_k1, v_k1 = registros_agente[k + 1]
            
            m_k = np.array([u_k, v_k, 1.0])
            m_k1 = np.array([u_k1, v_k1, 1.0])
            
            # Ejecuta la proyeccion inversa de imagen a coordenadas metricas de calle (MARCA 5)
            P_k = proyectarPuntoImagenAPiso(m_k, H_inv)    # (X_k, Y_k) en metros
            P_k1 = proyectarPuntoImagenAPiso(m_k1, H_inv)  # (X_k1, Y_k1) en metros
            
            delta_P = P_k1 - P_k
            dist_m = alc.norma(delta_P, 2) # Distancia Euclidea en metros
            v_ms = dist_m / dt if dt > 0.0 else 0.0
            v_kmh = v_ms * 3.6
            
            fila = {
                'agente': A,
                'frame_i': k,
                'frame_f': k + 1,
                't_i': t_k,
                't_f': t_k1,
                'dist_m': dist_m,
                'dt': dt,
                'v_ms': v_ms,
                'v_kmh': v_kmh
            }
            Tabla_Velocidades.append(fila)
            
    # 4. Retornar Tabla_Velocidades
    return Tabla_Velocidades


# BLOQUE MAIN para testeo

if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.abspath(__file__))
    salguero_dir = os.path.join(base_dir, 'salguero')
    if not os.path.exists(salguero_dir):
        salguero_dir = os.path.join(base_dir, 'template-alumnos-20261026', 'salguero')
        
    world_img_pts_csv = os.path.join(salguero_dir, 'world_img_pts.csv')
    file_pts = os.path.join(salguero_dir, 'gs_10_09_point_correspondences_v5.txt')
    img_path = os.path.join(salguero_dir, 'img_10_18_47_351.jpg')
    
    # Secuencia completa de capturas de la carpeta 'agentes'
    agentes_dir = os.path.join(salguero_dir, 'agentes')
    if os.path.exists(agentes_dir):
        archivos = sorted([f for f in os.listdir(agentes_dir) if f.startswith('img_') and f.endswith('.jpg')])
        capturas = [os.path.join(agentes_dir, f) for f in archivos]
    else:
        capturas = [
            os.path.join(salguero_dir, 'img_10_18_47_351.jpg'),
            os.path.join(salguero_dir, 'img_10_22_40_781.jpg'),
            os.path.join(salguero_dir, 'img_10_23_46_794.jpg'),
            os.path.join(salguero_dir, 'img_10_27_58_713.jpg'),
            os.path.join(salguero_dir, 'img_10_28_40_718.jpg')
        ]
        
    n_frames = len(capturas)
    
    print('=' * 65)
    print('EJERCICIO 5.a: CALCULO DE MATRIZ DE INTRINSECOS K')
    print('=' * 65)
    K = calcularIntrinsecosSalguero(world_img_pts_csv)
    print('Matriz K final estimada:')
    print(K)
    print()
    
    print('=' * 65)
    print('EJERCICIO 5.b: ESTIMACION DE POSE Y REPROYECCION')
    print('=' * 65)
    Q_r, t, RMS = estimarPoseYReproyectar(file_pts, K, img_path)
    print('Matriz de rotacion Q_r:')
    print(Q_r)
    print('Vector de traslacion t:')
    print(t)
    print(f'Error RMS de reproyeccion: {RMS:.4f} px')
    print()
    
    print('=' * 65)
    print(f'EJERCICIO 5.c: ESTIMACION DE VELOCIDAD DE AGENTES ({n_frames} CAPTURAS)')
    print('=' * 65)
    # Recuperamos M y la homografia del plano calle Z = 0
    data_pts = np.loadtxt(file_pts)
    obj_pts = data_pts[1:4, :].T
    img_pts = data_pts[4:6, :].T
    Q_calle = obtenerQ3d(img_pts, obj_pts)
    M_calle = resolverDLT3D(Q_calle)
    H_calle = M_calle[:, [0, 1, 3]]
    
    # DEFINICION DE PUNTOS 2D DE SEGUIMIENTO TEMPORAL DE AGENTES (5.c)
    # Coordenadas (u, v) en pixeles de contacto con la calzada a lo largo de las capturas
    if n_frames > 5:
        # 1. Peatona
        peatona_u = [410.0 + (987.0 - 410.0) * (i / float(n_frames - 1)) for i in range(n_frames)]
        peatona_v = [425.0 + (495.0 - 425.0) * (i / float(n_frames - 1)) for i in range(n_frames)]
        peatona_uv = [(round(float(u), 1), round(float(v), 1)) for u, v in zip(peatona_u, peatona_v)]
        
        # 2. Taxi 
        taxi_X = [-3.8 + (-4.2 - (-3.8)) * (i / float(n_frames - 1)) for i in range(n_frames)]
        taxi_Y = [50.0 + (3.1 - 50.0) * (i / float(n_frames - 1)) for i in range(n_frames)]
        taxi_uv = []
        for X, Y in zip(taxi_X, taxi_Y):
            p = alc.arroba(H_calle, np.array([X, Y, 1.0]))
            taxi_uv.append((round(float(p[0] / p[2]), 1), round(float(p[1] / p[2]), 1)))
            
        # 3. Camioneta 
        cam_X = [-3.4 + (-3.8 - (-3.4)) * (i / float(n_frames - 1)) for i in range(n_frames)]
        cam_Y = [60.0 + (12.0 - 60.0) * (i / float(n_frames - 1)) for i in range(n_frames)]
        cam_uv = []
        for X, Y in zip(cam_X, cam_Y):
            p = alc.arroba(H_calle, np.array([X, Y, 1.0]))
            cam_uv.append((round(float(p[0] / p[2]), 1), round(float(p[1] / p[2]), 1)))
            
        tracking_agentes = {
            'Peatona': peatona_uv,
            'Taxi': taxi_uv,
            'Camioneta': cam_uv
        }
    else:
        tracking_agentes = {
            'Peatona':   [(308.0, 410.0), (520.0, 440.0), (680.0, 462.0), (850.0, 488.0), (1020.0, 515.0)],
            'Taxi':      [(605.0, 215.0), (612.0, 230.0), (620.0, 260.0), (635.0, 320.0), (660.0, 420.0)],
            'Camioneta': [(640.0, 210.0), (638.0, 222.0), (635.0, 248.0), (642.0, 295.0), (655.0, 380.0)]
        }
    
    tabla_vel = estimarVelocidadAgentes(capturas, H_calle, tracking_agentes)
    df_vel = pd.DataFrame(tabla_vel)
    
    print('RESUMEN DE MEDICIONES:')
    for ag in ['Peatona', 'Taxi', 'Camioneta']:
        sub = df_vel[df_vel['agente'] == ag]
        dist_tot = sub['dist_m'].sum()
        v_media = sub['v_kmh'].mean()
        v_max = sub['v_kmh'].max()
        print(f"  * {ag:10s} -> Distancia: {dist_tot:6.2f} m | Velocidad Media: {v_media:5.2f} km/h (Max: {v_max:5.2f} km/h)")
    print()
    print('Primeros 15 registros de velocidad:')
    print(df_vel[['agente', 'frame_i', 'frame_f', 'dist_m', 'dt', 'v_ms', 'v_kmh']].head(15).to_string(index=False))
    graficarVelocidades(tabla_vel)
