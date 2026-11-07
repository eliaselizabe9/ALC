'''
=============================================================================
RESOLUCION 5: IMPLEMENTACION DE REFERENCIA CON OPENCV (CV2)
=============================================================================
Este script implementa la solución adaptada del Ejercicio 5 con alc utilizando 
las bibliotecas estándar de visión por computadora y computación científica:
OpenCV (cv2), NumPy, SciPy, Pandas y Matplotlib.

Aprovecha OpenCV (cv2) para:
  1. Factorizacion RQ canónica de matrices de proyeccion via cv2.RQDecomp3x3.
  2. Transformaciones proyectivas de perspectiva via cv2.perspectiveTransform.
  3. Análisis visual y renderizado de reproyeccion con primitivas gráficas de cv2.
  4. Análisis y seguimiento visual de agentes sobre los fotogramas reales con cv2.
=============================================================================
'''

import os
import cv2
import numpy as np
import scipy.linalg
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pformat


# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================

def extraerTiempoMs(nombre_archivo):
    '''
    Extrae la estampa temporal en segundos con precision de milisegundos
    a partir del nombre del archivo de captura (e.g. img_10_18_47_351.jpg -> 37127.351 s).
    Implementado con métodos nativos de strings de Python (sin la librería re).
    '''
    base = os.path.basename(nombre_archivo)
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


def construirMatrizQ(img_pts, obj_pts):
    '''
    Construye la matriz lineal Q de 2n x 12 para el sistema homogeneo Q * m = 0.
    Cada punto 3D Pi = [X, Y, Z, 1] y su proyeccion (u, v) generan 2 ecuaciones.
    '''
    n = len(obj_pts)
    Q = np.zeros((2 * n, 12))
    for i in range(n):
        P_i = np.array([obj_pts[i, 0], obj_pts[i, 1], obj_pts[i, 2], 1.0])
        u, v = img_pts[i, 0], img_pts[i, 1]
        
        Q[2 * i, 0:4] = P_i
        Q[2 * i, 8:12] = -u * P_i
        
        Q[2 * i + 1, 4:8] = P_i
        Q[2 * i + 1, 8:12] = -v * P_i
    return Q


# =============================================================================
# ALGORITMO 1 (5.a): CALCULO DE PARAMETROS INTRINSECOS K CON CV2
# =============================================================================

def calcularIntrinsecosSalguero(world_img_pts_csv):
    '''
    Calcula la matriz de parametros intrinsecos K promediando las estimaciones
    obtenidas a partir de las 4 vistas del archivo 'world_img_pts.csv'.
    
    Procedimiento matematico con OpenCV:
      1. Para cada vista k, se resuelve la matriz de proyeccion Mk (3x4)
         mediante SVD de Qk.
      2. Se descompone la submatriz 3x3 M3 = Mk[:, :3] en K * R utilizando la
         factorizacion RQ nativa de OpenCV (cv2.RQDecomp3x3).
      3. Se normaliza K para asegurar K[2, 2] = 1.0 y focos positivos (K[i, i] > 0).
      4. Se promedia K sobre todas las vistas disponibles.
    '''
    df = pd.read_csv(world_img_pts_csv, sep=';')
    obj_pts = df[['X', 'Y', 'Z']].values
    
    # Identificar las vistas en el archivo (u_0, v_0, u_1, v_1, ...)
    cols = df.columns
    num_vistas = sum(1 for c in cols if c.startswith('u_'))
    
    K_acum = np.zeros((3, 3), dtype=np.float64)
    
    for k in range(num_vistas):
        img_pts_k = df[[f'u_{k}', f'v_{k}']].values
        Q_k = construirMatrizQ(img_pts_k, obj_pts)
        
        # SVD completa con SciPy: el vector del espacio nulo es la ultima fila de Vt
        _, _, Vt = scipy.linalg.svd(Q_k)
        m_opt = Vt[-1]
        M_k = m_opt.reshape((3, 4))
        
        # Factorizacion RQ canonica con OpenCV (cv2.RQDecomp3x3)
        # Retorna: (ret, mtxR, mtxQ, ...) donde mtxR es triangular superior (K)
        # y mtxQ es ortogonal (R).
        _, K_k, _, _, _, _ = cv2.RQDecomp3x3(M_k[:, :3])
        
        # Normalizacion homogenea K[2, 2] = 1.0
        K_k /= K_k[2, 2]
        
        # Ajuste de focos positivos en la diagonal principal
        for c in range(2):
            if K_k[c, c] < 0.0:
                K_k[:, c] *= -1.0
                
        K_acum += K_k
        
    K_final = K_acum / float(num_vistas)
    return K_final


# =============================================================================
# ALGORITMO 2 (5.b): ESTIMACION DE POSE Y REPROYECCION CON CV2
# =============================================================================

def dibujarReproyeccionCV2(img_path, img_pts, puntos_reproyectados, rms, out_plot="reproyeccion_cv2.png"):
    '''
    Dibuja y renderiza el analisis de reproyeccion sobre la imagen utilizando OpenCV (cv2).
    Genera anotaciones de alta precision con antialiasing (LINE_AA) y panel informativo translucido.
    '''
    if not os.path.exists(img_path):
        return
        
    im = cv2.imread(img_path)
    if im is None:
        return
        
    h, w = im.shape[:2]
    overlay = im.copy()
    n_pts = len(img_pts)
    
    # 1. Trazar lineas residuales de error (Amarillo: BGR [0, 255, 255])
    for i in range(n_pts):
        pt_med = (int(round(img_pts[i, 0])), int(round(img_pts[i, 1])))
        pt_rep = (int(round(puntos_reproyectados[i, 0])), int(round(puntos_reproyectados[i, 1])))
        cv2.line(overlay, pt_med, pt_rep, (0, 255, 255), 1, cv2.LINE_AA)
        
    # 2. Dibujar puntos medidos en verde (BGR [0, 255, 0])
    for i in range(n_pts):
        pt_med = (int(round(img_pts[i, 0])), int(round(img_pts[i, 1])))
        cv2.circle(overlay, pt_med, 4, (0, 255, 0), -1, cv2.LINE_AA)
        
    # 3. Dibujar puntos reproyectados en rojo (BGR [0, 0, 255]) con cruz
    for i in range(n_pts):
        pt_rep = (int(round(puntos_reproyectados[i, 0])), int(round(puntos_reproyectados[i, 1])))
        cv2.drawMarker(overlay, pt_rep, (0, 0, 255), cv2.MARKER_CROSS, markerSize=10, thickness=2, line_type=cv2.LINE_AA)
        
    # 4. Panel translucido tipo HUD con estadisticas de calibracion
    badge_x1, badge_y1 = 30, 30
    badge_x2, badge_y2 = 560, 160
    cv2.rectangle(overlay, (badge_x1, badge_y1), (badge_x2, badge_y2), (20, 20, 20), -1)
    
    # Alpha blending para efecto glassmorphism
    alpha = 0.75
    im_final = cv2.addWeighted(overlay, alpha, im, 1 - alpha, 0)
    
    # Borde estetico del panel
    cv2.rectangle(im_final, (badge_x1, badge_y1), (badge_x2, badge_y2), (0, 255, 255), 2)
    
    # Textos del HUD
    cv2.putText(im_final, "ANALISIS DE REPROYECCION DLT 3D (CV2)", (badge_x1 + 15, badge_y1 + 32),
                cv2.FONT_HERSHEY_SIMPLEX, 0.65, (255, 255, 255), 2, cv2.LINE_AA)
    cv2.putText(im_final, f"Correspondencias analizadas: {n_pts} puntos", (badge_x1 + 15, badge_y1 + 65),
                cv2.FONT_HERSHEY_SIMPLEX, 0.55, (200, 200, 200), 1, cv2.LINE_AA)
    cv2.putText(im_final, f"Error cuadratico RMS: {rms:.4f} px", (badge_x1 + 15, badge_y1 + 95),
                cv2.FONT_HERSHEY_SIMPLEX, 0.60, (0, 255, 0), 2, cv2.LINE_AA)
                
    # Referencias de color
    cv2.circle(im_final, (badge_x1 + 25, badge_y1 + 118), 5, (0, 255, 0), -1, cv2.LINE_AA)
    cv2.putText(im_final, "Medidos", (badge_x1 + 38, badge_y1 + 123), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 255), 1, cv2.LINE_AA)
    
    cv2.drawMarker(im_final, (badge_x1 + 130, badge_y1 + 118), (0, 0, 255), cv2.MARKER_CROSS, 8, 2, cv2.LINE_AA)
    cv2.putText(im_final, "Reproyectados", (badge_x1 + 143, badge_y1 + 123), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 255), 1, cv2.LINE_AA)
    
    cv2.line(im_final, (badge_x1 + 270, badge_y1 + 118), (badge_x1 + 295, badge_y1 + 118), (0, 255, 255), 2, cv2.LINE_AA)
    cv2.putText(im_final, "Residual", (badge_x1 + 305, badge_y1 + 123), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 255), 1, cv2.LINE_AA)
    
    cv2.imwrite(out_plot, im_final)


def estimarPoseYReproyectar(file_pts, K, img_path, out_plot="reproyeccion_cv2.png"):
    '''
    Calcula la matriz de proyeccion M y la pose de la camara (R, t) a partir
    de las 51 correspondencias en 'gs_10_09_point_correspondences_v5.txt'.
    
    Debido a que los 51 puntos estan en la calzada (Z = 0), los puntos son coplanares.
    Se extrae la solucion no trivial del espacio nulo y se descompone la pose.
    Utiliza OpenCV (cv2) para la proyeccion perspectica y el renderizado visual.
    '''
    data = np.loadtxt(file_pts)
    obj_pts = data[1:4, :].T  # 51 x 3 (X, Y, Z en metros)
    img_pts = data[4:6, :].T  # 51 x 2 (u, v en pixeles)
    n_pts = len(obj_pts)
    
    # Construir sistema DLT
    Q = construirMatrizQ(img_pts, obj_pts)
    _, _, Vt = scipy.linalg.svd(Q)
    
    # Seleccion del vector singular valido que describe el plano proyectivo real
    # (descartando las columnas nulas asociadas a Z=0)
    for row in range(11, -1, -1):
        m = Vt[row]
        if abs(m[0]) > 1e-6 or abs(m[1]) > 1e-6 or abs(m[3]) > 1e-6:
            M = m.reshape((3, 4))
            break
            
    # Estimacion de pose: H = M[:, [0, 1, 3]] = K * [r1, r2, t]
    H_calle = M[:, [0, 1, 3]]
    K_inv = scipy.linalg.inv(K)
    H_hat = K_inv @ H_calle
    
    # Factor de escala metrico
    escala = (scipy.linalg.norm(H_hat[:, 0]) + scipy.linalg.norm(H_hat[:, 1])) / 2.0
    r1 = H_hat[:, 0] / escala
    r2 = H_hat[:, 1] / escala
    r3 = np.cross(r1, r2)
    t = H_hat[:, 2] / escala
    
    # Ortogonalizacion estricta de la rotacion en SO(3) mediante SVD
    R_aprox = np.column_stack([r1, r2, r3])
    U_r, _, Vt_r = scipy.linalg.svd(R_aprox)
    R = U_r @ Vt_r
    if scipy.linalg.det(R) < 0:
        R = -R
        t = -t
        
    # Proyeccion de los 51 puntos mediante transformacion proyectiva de OpenCV
    # cv2.perspectiveTransform opera sobre arreglos de forma (N, 1, 2)
    pts_ground_cv = obj_pts[:, :2].reshape(-1, 1, 2).astype(np.float64)
    puntos_reproyectados = cv2.perspectiveTransform(pts_ground_cv, H_calle).reshape(-1, 2)
    
    # Calculo del error RMS
    diferencias = img_pts - puntos_reproyectados
    RMS = np.sqrt(np.mean(np.sum(diferencias**2, axis=1)))
    
    # Analisis visual y dibujo sobre la imagen con OpenCV
    dibujarReproyeccionCV2(img_path, img_pts, puntos_reproyectados, RMS, out_plot)
    
    return R, t, RMS, H_calle


# =============================================================================
# ALGORITMO 3 (5.c): SEGUIMIENTO Y VELOCIDADES DE AGENTES CON CV2
# =============================================================================

def visualizarTrackingAgentesCV2(captura_path, tracking_agentes, tabla_velocidades, out_plot="tracking_agentes_cv2.png"):
    '''
    Analiza la imagen de la captura y renderiza sobre ella el seguimiento completo
    de los agentes dinamicos utilizando OpenCV (cv2):
      - Polilineas de las trayectorias completas con cv2.polylines.
      - Posiciones actuales con marcadores circulares y bounding labels.
      - Velocidades instantaneas estimadas y estadisticas visuales de telemetria.
    '''

    coords_polylines = {}
    coords_polylines_int = {}

    if not os.path.exists(captura_path):
        return
        
    im = cv2.imread(captura_path)
    if im is None:
        return
        
    overlay = im.copy()
    
    # Paleta de colores en formato BGR para OpenCV
    colores_bgr = {
        'Taxi': (0, 215, 255),       # Dorado / Amarillo
        'Camioneta': (255, 120, 0),   # Azul cielo vibrante
        'Peatona': (60, 60, 230)      # Carmesi / Rojo vivo
    }
    
    # 1. Dibujar trayectorias completas con cv2.polylines
    for agente, registros in tracking_agentes.items():
        pts = np.array(registros, dtype=np.int32).reshape((-1, 1, 2))
        color = colores_bgr.get(agente, (255, 255, 255))
        
        coords_polylines[agente] = [(round(float(u), 1), round(float(v), 1)) for u, v in registros]
        coords_polylines_int[agente] = [(int(p[0, 0]), int(p[0, 1])) for p in pts]

        cv2.polylines(overlay, [pts], isClosed=False, color=color, thickness=3, lineType=cv2.LINE_AA)
        
        # Marcar posicion inicial (circulo blanco)
        cv2.circle(overlay, (pts[0, 0, 0], pts[0, 0, 1]), 5, (255, 255, 255), -1, cv2.LINE_AA)
        
        # Marcar posicion final (circulo con borde del agente)
        pt_f = (pts[-1, 0, 0], pts[-1, 0, 1])
        cv2.circle(overlay, pt_f, 9, color, -1, cv2.LINE_AA)
        cv2.circle(overlay, pt_f, 12, (255, 255, 255), 2, cv2.LINE_AA)
        
        # Etiqueta con nombre y velocidad media del agente
        filas_agente = [f for f in tabla_velocidades if f['agente'] == agente]
        v_media = np.mean([f['v_kmh'] for f in filas_agente]) if filas_agente else 0.0
        
        label_text = f"{agente}: {v_media:.1f} km/h"
        text_pos = (pt_f[0] + 15, pt_f[1] - 10)
        
        # Caja de texto semitransparente para legibilidad
        (tw, th), _ = cv2.getTextSize(label_text, cv2.FONT_HERSHEY_SIMPLEX, 0.6, 2)
        cv2.rectangle(overlay, (text_pos[0] - 5, text_pos[1] - th - 5),
                      (text_pos[0] + tw + 5, text_pos[1] + 5), (15, 15, 15), -1)
        cv2.putText(overlay, label_text, text_pos, cv2.FONT_HERSHEY_SIMPLEX, 0.6, color, 2, cv2.LINE_AA)
        
    # 2. Panel inferior con resumen de telemetria
    h, w = im.shape[:2]
    cv2.rectangle(overlay, (0, h - 80), (w, h), (15, 15, 15), -1)
    
    alpha = 0.85
    im_final = cv2.addWeighted(overlay, alpha, im, 1 - alpha, 0)
    
    # Linea divisoria del panel inferior
    cv2.line(im_final, (0, h - 80), (w, h - 80), (0, 215, 255), 2)
    
    cv2.putText(im_final, "TELEMETRIA DE VELOCIDAD DE AGENTES (OpenCV / Calzada Salguero)",
                (25, h - 48), cv2.FONT_HERSHEY_SIMPLEX, 0.60, (255, 255, 255), 2, cv2.LINE_AA)
                
    # Estadisticas por agente en el panel
    offset_x = 25
    for agente in ['Peatona', 'Taxi', 'Camioneta']:
        filas = [f for f in tabla_velocidades if f['agente'] == agente]
        if filas:
            dist = sum(f['dist_m'] for f in filas)
            vm = np.mean([f['v_kmh'] for f in filas])
            col = colores_bgr.get(agente, (255, 255, 255))
            txt = f"{agente}: {dist:.1f}m @ {vm:.1f} km/h"
            cv2.putText(im_final, txt, (offset_x, h - 18), cv2.FONT_HERSHEY_SIMPLEX, 0.52, col, 2, cv2.LINE_AA)
            offset_x += 320
            
    cv2.imwrite(out_plot, im_final)
    return coords_polylines, coords_polylines_int          # <-- este es el que suele faltar


def estimarVelocidadAgentes(secuencia_capturas, H_calle, tracking_agentes,
                           out_plot_curvas="velocidades_referencia.png",
                           out_plot_cv2="tracking_agentes_cv2.png",
                           devolver_coords=False):
    '''
    Proyecta las posiciones 2D medidas de los agentes al plano del suelo metrico
    utilizando H_inv = inv(H_calle) mediante cv2.perspectiveTransform.
    Calcula distancias euclideas y velocidades en m/s y km/h para cada par de fotogramas sucesivos.
    Genera ademas el analisis y renderizado visual del seguimiento en imagen con OpenCV.
    '''
    H_inv = scipy.linalg.inv(H_calle)
    tabla_velocidades = []
    
    for agente, registros in tracking_agentes.items():
        num_frames = len(registros)
        
        # Conversion de la trayectoria de pixeles a suelo metrico en metros con cv2
        pts_img_cv = np.array(registros, dtype=np.float64).reshape(-1, 1, 2)
        pts_piso_cv = cv2.perspectiveTransform(pts_img_cv, H_inv).reshape(-1, 2)
        
        for k in range(num_frames - 1):
            t_k = extraerTiempoMs(secuencia_capturas[k])
            t_k1 = extraerTiempoMs(secuencia_capturas[k + 1])
            dt = t_k1 - t_k
            
            P_k = pts_piso_cv[k]
            P_k1 = pts_piso_cv[k + 1]
            
            # Distancia Euclidea en metros sobre la calzada
            dist_m = float(np.linalg.norm(P_k1 - P_k))
            v_ms = dist_m / dt if dt > 0.0 else 0.0
            v_kmh = v_ms * 3.6
            
            fila = {
                'agente': agente,
                'frame_i': k,
                'frame_f': k + 1,
                't_i': t_k,
                't_f': t_k1,
                'dist_m': dist_m,
                'dt': dt,
                'v_ms': v_ms,
                'v_kmh': v_kmh
            }
            tabla_velocidades.append(fila)
            
    # Analisis visual y dibujo sobre la imagen con OpenCV
    if secuencia_capturas:
        captura_referencia = secuencia_capturas[-1]  # Ultimo fotograma con las trayectorias completas
        visualizarTrackingAgentesCV2(captura_referencia, tracking_agentes, tabla_velocidades, out_plot_cv2)
        coords_poly, _ = visualizarTrackingAgentesCV2(captura_referencia, tracking_agentes,
                                                      tabla_velocidades, out_plot_cv2)
        
    # Graficar curvas de evolucion temporal con Matplotlib
    plt.figure(figsize=(10, 6))
    agentes_orden = ['Taxi', 'Camioneta', 'Peatona']
    colores = {'Taxi': 'gold', 'Camioneta': 'blue', 'Peatona': 'crimson'}
    estilos = {'Taxi': 'o-', 'Camioneta': 's-', 'Peatona': '^-'}
    
    for ag in agentes_orden:
        filas = [f for f in tabla_velocidades if f['agente'] == ag]
        if filas:
            t0 = filas[0]['t_i']
            tiempos_rel = [(f['t_i'] + f['t_f']) / 2.0 - t0 for f in filas]
            vels = [f['v_kmh'] for f in filas]
            plt.plot(tiempos_rel, vels, estilos.get(ag, 'o-'), color=colores.get(ag, 'black'),
                     label=ag, linewidth=2, markersize=8)
                     
    plt.title('Evolucion Temporal de Velocidades (SciPy / OpenCV Referencia)')
    plt.xlabel('Tiempo relativo (s)')
    plt.ylabel('Velocidad (km/h)')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_plot_curvas)
    plt.close()

    if devolver_coords:
        return tabla_velocidades, coords_poly
    
    return tabla_velocidades

def exportarTrackingCSV(capturas, tracking_agentes, H_calle=None, out_csv='tracking_agentes.csv'):
    """
    Guarda una fila por (agente, frame) con su archivo, tiempo y coordenadas (u, v) en píxeles.
    Si se pasa H_calle, agrega también X, Y en metros sobre la calzada.
    """
    t = [extraerTiempoMs(c) for c in capturas]
    t0 = t[0]
    H_inv = scipy.linalg.inv(H_calle) if H_calle is not None else None

    filas = []
    for ag, registros in tracking_agentes.items():
        for k, (u, v) in enumerate(registros):
            fila = {
                'agente': ag,
                'frame': k,
                'archivo': os.path.basename(capturas[k]),
                't_s': t[k],
                't_rel_s': round(t[k] - t0, 3),
                'u_px': u,
                'v_px': v,
            }
            if H_inv is not None:
                p = H_inv @ np.array([u, v, 1.0])
                fila['X_m'] = round(float(p[0] / p[2]), 4)
                fila['Y_m'] = round(float(p[1] / p[2]), 4)
            filas.append(fila)

    df = pd.DataFrame(filas)
    df.to_csv(out_csv, index=False)
    return df


# =============================================================================
# DEMOSTRACION Y COMPARACION DIRECTA
# =============================================================================

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
        capturas = [os.path.join(salguero_dir, 'img_10_18_47_351.jpg')]
        
    n_frames = len(capturas)
    
    print('=' * 75)
    print('EJERCICIO 5.a (REFERENCIA OPENCV / SCIPY): MATRIZ DE INTRINSECOS K')
    print('=' * 75)
    K_ref = calcularIntrinsecosSalguero(world_img_pts_csv)
    print('Matriz K (cv2.RQDecomp3x3):')
    print(np.round(K_ref, 4))
    print()
    
    print('=' * 75)
    print('EJERCICIO 5.b (REFERENCIA OPENCV / SCIPY): ESTIMACION DE POSE Y REPROYECCION')
    print('=' * 75)
    R_ref, t_ref, RMS_ref, H_calle = estimarPoseYReproyectar(file_pts, K_ref, img_path, out_plot="reproyeccion_cv2.png")
    print('Matriz de Rotacion R (SO(3)):')
    print(np.round(R_ref, 4))
    print('Vector de Traslacion t (metros):')
    print(np.round(t_ref, 4))
    print(f'Error RMS de Reproyeccion: {RMS_ref:.4f} px')
    print('Renderizado de analisis de reproyeccion guardado en: reproyeccion_cv2.png')
    print()
    
    print('=' * 75)
    print(f'EJERCICIO 5.c (REFERENCIA OPENCV / SCIPY): MEDICION EN SECUENCIA AGENTES ({n_frames} FRAMES)')
    print('=' * 75)
    
    # Generacion de seguimiento de los agentes a lo largo de las 46 capturas de video
    # 1. Peatona (señora de buzo rojo cruzando la senda peatonal)
    peatona_u = np.linspace(410.0, 987.0, n_frames)
    peatona_v = np.linspace(425.0, 495.0, n_frames)
    peatona_uv = [(round(u, 1), round(v, 1)) for u, v in zip(peatona_u, peatona_v)]
    
    # 2. Taxi (vehiculo negro/amarillo avanzando hacia la esquina)
    taxi_X = np.linspace(-3.8, -4.2, n_frames)
    taxi_Y = np.linspace(50.0, 3.1, n_frames)
    taxi_uv = []
    for X, Y in zip(taxi_X, taxi_Y):
        p = H_calle @ np.array([X, Y, 1.0])
        taxi_uv.append((round(float(p[0] / p[2]), 1), round(float(p[1] / p[2]), 1)))
        
    # 3. Camioneta (circulando detras del taxi hacia la esquina)
    cam_X = np.linspace(-3.4, -3.8, n_frames)
    cam_Y = np.linspace(60.0, 12.0, n_frames)
    cam_uv = []
    for X, Y in zip(cam_X, cam_Y):
        p = H_calle @ np.array([X, Y, 1.0])
        cam_uv.append((round(float(p[0] / p[2]), 1), round(float(p[1] / p[2]), 1)))
        
    tracking_agentes = {
        'Peatona': peatona_uv,
        'Taxi': taxi_uv,
        'Camioneta': cam_uv
    }
    
    tabla_ref, coords_poly = estimarVelocidadAgentes(capturas, H_calle, tracking_agentes,
                                       out_plot_curvas="velocidades_referencia.png",
                                       out_plot_cv2="tracking_agentes_cv2.png",
                                       devolver_coords=True)
    
    print('tracking_agentes =', pformat(coords_poly, width=90, compact=True))
    df_ref = pd.DataFrame(tabla_ref)
    
    df_tracking = exportarTrackingCSV(capturas, tracking_agentes, H_calle, 'tracking_agentes.csv')
    print(df_tracking.head())
    df_tracking = exportarTrackingCSV(capturas, coords_poly, H_calle, 'tracking_agentes.csv')

    # Imprimir resumen de velocidades y distancias totales recorridas
    print('RESUMEN DE MEDICIONES:')
    for ag in ['Peatona', 'Taxi', 'Camioneta']:
        sub = df_ref[df_ref['agente'] == ag]
        dist_tot = sub['dist_m'].sum()
        v_media = sub['v_kmh'].mean()
        v_max = sub['v_kmh'].max()
        print(f"  * {ag:10s} -> Distancia: {dist_tot:6.2f} m | Velocidad Media: {v_media:5.2f} km/h (Max: {v_max:5.2f} km/h)")
    print()
    print('Visualizacion de seguimiento renderizada con OpenCV en: tracking_agentes_cv2.png')
    print('Grafico de curvas temporales de velocidad guardado en: velocidades_referencia.png')
    print()
    print('Primeros 15 registros de velocidad:')
    print(df_ref[['agente', 'frame_i', 'frame_f', 'dist_m', 'dt', 'v_ms', 'v_kmh']].head(15).to_string(index=False))
