# Notación y Convenciones Fundamentales
Consideramos el sistema de ecuaciones lineales $Ax = b$ con $A \in \mathbb{R}^{n \times n}$ invertible y $b \in \mathbb{R}^n$. A lo largo de todo el documento se mantiene la siguiente notación uniforme:

|**Símbolo**|**Significado**|**Definición / Propiedad**|
|---|---|---|
|$x^*$|Solución exacta del sistema|$x^* = A^{-1}b$|
|$x^{(k)}$|Vector aproximación en la iteración $k$|$x^{(0)}$ es el vector inicial arbitrario|
|$e^{(k)}$|Vector de error en el paso $k$|$e^{(k)} = x^{(k)} - x^*$ (o $x^* - x^{(k)}$ según convención; aquí usamos $e^{(k)} = x^{(k)} - x^*$)|
|$r^{(k)}$|Vector residuo en el paso $k$|$r^{(k)} = b - Ax^{(k)}$|
|$\rho(G)$|Radio espectral de la matriz $G$|$\rho(G) = \max \{|\lambda| : \lambda \in \sigma(G)\}$|
|$\kappa_2(A)$|Número de condición espectral|$|A|_2 |A^{-1}|_2 = \sigma_{\max}/\sigma_{\min}$ (Si $A$ es SPD, $\kappa_2(A) = \lambda_{\max}/\lambda_{\min}$)|
|$\langle u, v \rangle$|Producto interno euclídeo canónico|$\langle u, v \rangle = u^T v$|
|$|v|_A$|Norma de energía (o $A$-norma, para $A$ SPD)|$|v|_A = \sqrt{v^T A v} = \sqrt{\langle v, v \rangle_A}$|

**Ecuación fundamental del defecto (Relación Error–Residuo):** Dado que $Ax^* = b$, se verifica en toda iteración que:
$$r^{(k)} = b - Ax^{(k)} = Ax^* - Ax^{(k)} = -A(x^{(k)} - x^*) = -Ae^{(k)} \iff Ae^{(k)} = -r^{(k)}$$
Mientras que el error $e^{(k)}$ es inaccesible en la práctica por desconocerse $x^*$, el residuo $r^{(k)}$ es directamente computable y mide cuánto falla la aproximación actual en satisfacer la ecuación.

# 1. Marco Teórico de Métodos Iterativos Clásicos (Kincaid 4.6 & Strang II.1)

A diferencia de los métodos directos (como la eliminación gaussiana $PA = LU$ o la factorización de Cholesky $A = LL^T$, que exigen $\mathcal{O}(n^3)$ operaciones y sufren el fenómeno de _fill-in_ o relleno de ceros en matrices ralas), los métodos iterativos generan una sucesión $\{x^{(k)}\}_{k=0}^\infty$ a partir de operaciones matriz-vector de costo $\mathcal{O}(\text{nnz}(A))$, preservando intacta la estructura de dispersión de $A$.

## 1.1. Formulación Matricial General mediante Partición (Splitting $Q$)

Siguiendo la formulación de Kincaid & Cheney (Sección 4.6) y Strang (Sección II.1), introducimos una **matriz de partición** (o _splitting matrix_) $Q \in \mathbb{R}^{n \times n}$ no singular, elegida de modo que los sistemas lineales de la forma $Qz = d$ sean computacionalmente triviales de resolver (por ejemplo, $Q$ diagonal o triangular).
Descomponiendo $A = Q - (Q - A)$ (en la notación de Strang, $A = S - T$ con $S = Q$ y $T = Q - A$), el sistema original $Ax = b$ se reescribe de manera equivalente como un problema de punto fijo:
$$Ax = b \iff Qx = (Q - A)x + b$$
Esto motiva el esquema iterativo estacionario lineal de primer orden:
$$Q x^{(k)} = (Q - A) x^{(k-1)} + b, \qquad k\geq1$$
Multiplicando a izquierda por $Q^{-1}$, obtenemos la forma explícita de la evolución iterativa:
$$x^{(k)} = (I - Q^{-1}A)x^{(k-1)} + Q^{-1}b = G x^{(k-1)} + c$$
donde definimos:
- **Matriz de iteración:** $G = I - Q^{-1}A = Q^{-1}(Q - A)$.
- **Vector de actualización constante:** $c = Q^{-1}b$.
 Notemos además que la corrección en cada paso responde directamente al residuo precondicionado por $Q^{-1}$:
$$x^{(k)} = x^{(k-1)} + Q^{-1}(b - Ax^{(k-1)}) = x^{(k-1)} + Q^{-1}r^{(k-1)}$$
## 1.2. Análisis de Convergencia Global: Condición Necesaria y Suficiente

Como la solución exacta $x^*$ satisface la ecuación de consistencia $Qx^* = (Q - A)x^* + b$ (es decir, $x^* = Gx^* + c$), restando ambas expresiones deducimos la ley de propagación del error:
$$Q(x^{(k)} - x^*) = (Q - A)(x^{(k-1)} - x^*) \implies e^{(k)} = G e^{(k-1)} = G^2 e^{(k-2)} = \dots = G^k e^{(0)}$$
>[!note] **Teorema de Convergencia Global de Métodos Estacionarios — Kincaid & Cheney 4.6:**
> La sucesión generada por $x^{(k)} = Gx^{(k-1)} + c$ converge a la única solución $x^* = A^{-1}b$ para **cualquier** vector inicial $x^{(0)} \in \mathbb{R}^n$ si y solo si el radio espectral de la matriz de iteración es estrictamente menor que la unidad:
> $$\rho(G) < 1$$
### Demostración
La condición $e^{(k)} = G^k e^{(0)} \to 0$ para todo $e^{(0)} \in \mathbb{R}^n$ es equivalente a afirmar que $\lim_{k \to \infty} G^k = 0$ (la matriz nula).
1. **Necesidad ($\implies$):** Supongamos que $\lim_{k \to \infty} e^{(k)} = 0$ para todo $e^{(0)} \in \mathbb{R}^n$. Sea $\lambda$ un autovalor cualquiera de $G$ y $v \neq 0$ su autovector asociado ($Gv = \lambda v$). Eligiendo el dato inicial tal que $e^{(0)} = v$, tenemos que $e^{(k)} = G^k v = \lambda^k v$. Como $v \neq 0$ y $\lim_{k \to \infty} \lambda^k v = 0$, necesariamente $|\lambda|^k \to 0$, lo cual exige $|\lambda| < 1$. Dado que esto vale para todo autovalor $\lambda \in \sigma(G)$, concluimos que $\rho(G) = \max_{\lambda} |\lambda| < 1$.
2. **Suficiencia ($\Longleftarrow$):** Supongamos que $\rho(G) < 1$. Por el Lema de Schur / Forma Canónica de Jordan (o el teorema de normas adaptadas de Kincaid, Sección 4.6), para todo $\epsilon > 0$ existe una norma vectorial $\|\cdot\|_\epsilon$ en $\mathbb{R}^n$ cuya norma matricial subordinada (inducida) verifica:$$|G|_\epsilon \le \rho(G) + \epsilon$$ Como $\rho(G) < 1$, podemos elegir $\epsilon = \frac{1 - \rho(G)}{2} > 0$, de modo que $\theta := \rho(G) + \epsilon = \frac{1 + \rho(G)}{2} < 1$. Aplicando la propiedad submultiplicativa de la norma inducida:$$|e^{(k)}|_\epsilon = |G^k e^{(0)}|_\epsilon \le |G|_\epsilon^k |e^{(0)}|_\epsilon \le \theta^k |e^{(0)}|_\epsilon$$Como $0 \le \theta < 1$, al tomar $k \to \infty$ resulta $\theta^k \to 0$, por lo que $\|e^{(k)}\|_\epsilon \to 0$. Por la equivalencia topológica de todas las normas en espacios de dimensión finita $\mathbb{R}^n$, $e^{(k)} \to 0$ en cualquier norma, luego $x^{(k)} \to x^*$ para todo $x^{(0)}$. $\blacksquare$
 **Tasa asintótica de convergencia:** La cantidad $R_\infty(G) = -\ln \rho(G)$ se denomina _velocidad asintótica de convergencia_. El número de iteraciones $k$ necesario para reducir la norma del error en un factor $10^{-m}$ se estima como $k \approx \frac{m \ln 10}{-\ln \rho(G)}$.
## 1.3. Desarrollo Detallado de los Métodos Clásicos
Descomponemos aditivamente la matriz $A$ en tres bloques estructurales:
$$A = L + D + U$$
donde $D = \text{diag}(a_{11}, \dots, a_{nn})$ es la diagonal principal (suponemos $a_{ii} \neq 0$), $L$ es la parte estrictamente triangular inferior ($l_{ij} = a_{ij}$ si $i > j$, $0$ en otro caso) y $U$ es la parte estrictamente triangular superior ($u_{ij} = a_{ij}$ si $i < j$, $0$ en otro caso). (Nota: En Kincaid 4.6 a veces se escala por la diagonal o se usan signos negativos $A = D - \tilde{L} - \tilde{U}$; aquí mantenemos la convención limpia $A = L + D + U$ del apunte base).
### 1. Método de Jacobi (Desplazamientos Simultáneos)
- **Elección de partición:** $Q_J = D$ y $Q_J - A = -(L + U)$.
- **Matriz de iteración:** $G_J = I - D^{-1}A = -D^{-1}(L + U)$.
- **Fórmula por componentes:** Se despeja la $i$-ésima incógnita utilizando exclusivamente las componentes del paso anterior $k-1$:$$x_i^{(k)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1, j \neq i}^n a_{ij} x_j^{(k-1)} \right), \qquad i = 1, \dots, n$$
- **Características computacionales:** Las $n$ ecuaciones son desacopladas dentro del mismo paso $k$. Es masivamente paralelizable (ideal para arquitecturas GPU/SIMD), pero exige almacenar en memoria dos vectores de estado ($x^{(k)}$ y $x^{(k-1)}$).
### 2. Método de Gauss-Seidel (Desplazamientos Sucesivos)
- **Elección de partición:** $Q_{GS} = D + L$ (triangular inferior) y $Q_{GS} - A = -U$. 
- **Matriz de iteración:** $G_{GS} = I - (D + L)^{-1}A = -(D + L)^{-1}U$. 
- **Fórmula por componentes:** Utiliza las componentes $x_1^{(k)}, \dots, x_{i-1}^{(k)}$ tan pronto como son calculadas en el mismo barrido: 
$$x_i^{(k)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i+1}^n a_{ij} x_j^{(k-1)} \right), \qquad i = 1, \dots, n$$
- **Características computacionales:** Sobrescribe el vector $x$ _in-place_ (requiere solo 1 vector de almacenamiento de estado), pero introduce dependencia de datos secuencial (sustitución hacia adelante), lo que dificulta su paralelización directa (salvo mediante técnicas de reordenamiento _Red-Black_ o coloreo de grafos).
### 3. Método SOR (Successive Over-Relaxation — Sobrerrelajación Sucesiva)
Para acelerar la convergencia de Gauss-Seidel, introducimos un parámetro de relajación $\omega \in \mathbb{R}$ que pondera el estado previo $x_i^{(k-1)}$ con la predicción pura de Gauss-Seidel $\tilde{x}_i^{(k)}$:
$$\begin{aligned}
	x_i^{(k)} &= (1 - \omega)x_i^{(k-1)} + \omega \tilde{x}_i^{(k)} \\
	&= x_i^{(k-1)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - a_{ii}x_i^{(k-1)} - \sum_{j=i+1}^n a_{ij} x_j^{(k-1)} \right)
\end{aligned}$$
Multiplicando por $a_{ii}/\omega$ y agrupando los términos en el paso $(k)$ a la izquierda y en el paso $(k-1)$ a la derecha:
$$\frac{1}{\omega} D x^{(k)} + L x^{(k)} = \left( \frac{1 - \omega}{\omega} D - U \right) x^{(k-1)} + b$$
- **Elección de partición:** $Q_{\text{SOR}} = \frac{1}{\omega}D + L = \frac{1}{\omega}(D + \omega L)$.
- **Matriz de iteración:** $G_\omega = (D + \omega L)^{-1} \big( (1 - \omega)D - \omega U \big)$.
- **Casos particulares:** Si $\omega = 1$, recuperamos Gauss-Seidel. Si $0 < \omega < 1$, se denomina _subrrelajación_ (útil para estabilizar sistemas no lineales). Si $\omega > 1$, se denomina _sobrerrelajación_.
#### Resultados fundamentales sobre SOR (Kincaid 4.6)
1. **Teorema de Kahan (Condición necesaria sobre $\omega$):** Para cualquier matriz $A$, el radio espectral de $G_\omega$ verifica $\rho(G_\omega) \ge |\omega - 1|$. En consecuencia, **SOR solo puede converger si $0 < \omega < 2$**.
	_Demostración:_ Los autovalores de $G_\omega$ tienen como producto $\det(G_\omega) = \det((D+\omega L)^{-1}) \det((1-\omega)D - \omega U) = (\det D)^{-1} (1-\omega)^n \det D = (1-\omega)^n$. Luego, $\rho(G_\omega)^n \ge |\det(G_\omega)| = |1-\omega|^n \implies \rho(G_\omega) \ge |\omega - 1|$. Para que $\rho(G_\omega) < 1$, es indispensable que $|\omega - 1| < 1 \iff 0 < \omega < 2$. $\blacksquare$
2. **Teorema de Ostrowski-Reich:** Si $A$ es simétrica y definida positiva (SPD), entonces SOR converge para todo $x^{(0)}$ si y solo si $0 < \omega < 2$ (en particular, Gauss-Seidel con $\omega=1$ siempre converge para matrices SPD).
3. **Parámetro óptimo de Young:** Si $A$ es consistentemente ordenada (como las matrices tridiagonales por bloques provenientes del laplaciano discretizado), los autovalores $\lambda$ de $G_\omega$ y $\mu$ de $G_J$ se relacionan por $(\lambda + \omega - 1)^2 = \lambda \omega^2 \mu^2$. El valor óptimo que minimiza $\rho(G_\omega)$ es:$$\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho(G_J)^2}}, \qquad \text{con } \rho(G_{\omega_{\text{opt}}}) = \omega_{\text{opt}} - 1$$ Para la ecuación de Poisson en una malla $N \times N$, donde $\rho(G_J) \approx 1 - \frac{\pi^2}{2N^2}$,el radio espectral pasa de $1 - \mathcal{O}(N^{-2})$ en Gauss-Seidel a $1 - \mathcal{O}(N^{-1})$ en SOR óptimo, reduciendo las iteraciones en un orden de magnitud completo.
### 4. Método SSOR (Symmetric Successive Over-Relaxation)
Cuando $A$ es simétrica ($U = L^T$), la matriz de partición de SOR ($Q_{\text{SOR}} = \frac{1}{\omega}D + L$) **no es simétrica**, lo que impide usar SOR como precondicionador simétrico o aplicarle aceleración polinomial directa de Chebyshev (que requiere autovalores reales). SSOR corrige esto ejecutando cada iteración en dos medios pasos: un barrido SOR hacia adelante ($i = 1 \dots n$) seguido de un barrido SOR hacia atrás ($i = n \dots 1$):
1. **Barrido hacia adelante:** $\frac{1}{\omega}(D + \omega L) x^{(k - 1/2)} = \frac{1}{\omega}\big((1 - \omega)D - \omega U\big) x^{(k-1)} + b$
2. **Barrido hacia atrás:** $\frac{1}{\omega}(D + \omega U) x^{(k)} = \frac{1}{\omega}\big((1 - \omega)D - \omega L\big) x^{(k - 1/2)} + b$
 Eliminando el vector intermedio $x^{(k - 1/2)}$, obtenemos la **matriz de partición de SSOR**:
$$Q_{\text{SSOR}} = \frac{1}{\omega(2 - \omega)} (D + \omega L) D^{-1} (D + \omega U)$$
Si $A$ es SPD ($U = L^T$ y $D$ positiva) y $0 < \omega < 2$, entonces $Q_{\text{SSOR}}$ es **simétrica y definida positiva**. Esto hace que la matriz de iteración $G_{\text{SSOR}} = I - Q_{\text{SSOR}}^{-1}A$ sea semejante a $I - Q_{\text{SSOR}}^{-1/2} A Q_{\text{SSOR}}^{-1/2}$ (simétrica), garantizando que todos los autovalores de $G_{\text{SSOR}}$ sean reales y positivos en $(0, 1)$. Por ello, $Q_{\text{SSOR}}$ es un excelente precondicionador simétrico (PCG) y es directamente compatible con la aceleración de Chebyshev.
## 1.4. Técnicas de Aceleración de la Convergencia (Kincaid 4.6)
Supongamos que disponemos de un método iterativo base $x^{(k)} = Gx^{(k-1)} + c$ donde los autovalores de $G$ son reales y están contenidos en un intervalo conocido $[\alpha, \beta]$ con $-1 < \alpha < \beta < 1$(o bien en $[a, b]$ incluso si $\ge 1$, siempre que $1 \notin [a, b]$).
### A. Extrapolación Analítica (Enriquecimiento en un paso)
Introducimos un parámetro $\gamma \neq 0$ para formar una combinación afín entre el estado previo y el paso propuesto por el método base:
$$x^{(k)} = \gamma \big( G x^{(k-1)} + c \big) + (1 - \gamma) x^{(k-1)} = G_\gamma x^{(k-1)} + \gamma c, \qquad \text{donde } G_\gamma = \gamma G + (1 - \gamma)I$$
Por el teorema de mapeo espectral, si $\lambda \in [\alpha, \beta]$ es un autovalor de $G$, el autovalor correspondiente de $G_\gamma$ es $\mu(\lambda) = \gamma \lambda + 1 - \gamma$. Queremos elegir $\gamma$ para minimizar $\rho(G_\gamma) = \max_{\lambda \in [\alpha, \beta]} |\gamma \lambda + 1 - \gamma|$. El equilibrio minimax óptimo ocurre cuando los extremos del intervalo mapeado son simétricos respecto del origen:
$$\gamma \alpha + 1 - \gamma = -(\gamma \beta + 1 - \gamma) \implies \gamma_{\text{opt}} = \frac{2}{2 - (\alpha + \beta)}$$
Con esta elección, el nuevo radio espectral se reduce a:
$$\rho(G_{\gamma_{\text{opt}}}) = \frac{\beta - \alpha}{2 - (\alpha + \beta)} < \max(|\alpha|, |\beta|)$$
### B. Aceleración Polinomial de Chebyshev (Semi-iteración de Golub-Varga)
En lugar de usar solo el último vector $x^{(k-1)}$, consideramos una combinación lineal de todas las iteraciones previas generadas por el método base $x^{(0)}, x^{(1)}, \dots, x^{(k)}$:
$$u^{(k)} = \sum_{j=0}^k a_{kj} x^{(j)}, \qquad \text{con la condición de consistencia } \sum_{j=0}^k a_{kj} = 1$$
Si $x^{(j)} = x^*$, entonces $u^{(k)} = x^*$. Para el error en $u^{(k)}$ tenemos:
$$u^{(k)} - x^* = \sum_{j=0}^k a_{kj} (x^{(j)} - x^*) = \left( \sum_{j=0}^k a_{kj} G^j \right) e^{(0)} = P_k(G) e^{(0)}$$
donde $P_k(t) = \sum_{j=0}^k a_{kj} t^j$ es un polinomio de grado $\le k$ normalizado tal que $P_k(1) = 1$. Si $G$ es simétrica (o semejante a una simétrica, como en Jacobi o SSOR sobre $A$ SPD), en cualquier norma adaptada tenemos:
$$|u^{(k)} - x^*| \le |P_k(G)| |e^{(0)}| = \left( \max_{\lambda \in [\alpha, \beta]} |P_k(\lambda)| \right) |e^{(0)}|$$
El problema de diseño óptimo consiste en hallar el polinomio $P_k \in \Pi_k$ con $P_k(1) = 1$ que minimice $\max_{t \in [\alpha, \beta]} |P_k(t)|$. Por la teoría clásica de aproximación (Kincaid, Sección 4.6 y 6.9), la solución única se obtiene reescalando el **polinomio de Chebyshev de primera especie** $T_k(z) = \cos(k \arccos z)$ desde $[-1, 1]$ al intervalo $[\alpha, \beta]$ mediante la transformación lineal $z(t) = \frac{2t - (\beta + \alpha)}{\beta - \alpha}$:
$$P_k(t) = \frac{T_k\left( \frac{2t - \beta - \alpha}{\beta - \alpha} \right)}{T_k\left( \frac{2 - \beta - \alpha}{\beta - \alpha} \right)} = \frac{T_k(z(t))}{T_k(z(1))}$$
donde $z(1) = \frac{2 - \beta - \alpha}{\beta - \alpha} > 1$. Dado que $|T_k(z)| \le 1$ para todo $z \in [-1, 1]$, el factor de atenuación del error está acotado por $\frac{1}{T_k(z(1))}$.
#### Recurrencia de tres términos libre de almacenamiento histórico:
Almacenar todos los vectores $x^{(0)}, \dots, x^{(k)}$ sería prohibitivo en memoria. Explotando la relación de recurrencia de los polinomios de Chebyshev, $T_{k+1}(z) = 2z T_k(z) - T_{k-1}(z)$ (con $T_0(z)=1, T_1(z)=z$), podemos calcular $u^{(k+1)}$ directamente a partir de $u^{(k)}$ y $u^{(k-1)}$. Definiendo $\rho_k = \frac{T_{k-1}(z(1))}{T_k(z(1))}$ y el parámetro de extrapolación $\gamma = \frac{2}{2 - \alpha - \beta}$, la recurrencia de aceleración de Chebyshev toma la forma:
1. Inicialización: $u^{(0)} = x^{(0)}$, y primer paso extrapolado: $u^{(1)} = \gamma(G u^{(0)} + c) + (1 - \gamma)u^{(0)}$.
 2. Para $k \ge 1$, los pesos dinámicos $\omega_{k+1}$ se actualizan por:
$$\omega_1 = 1, \quad \omega_2 = \frac{2 z(1)^2}{2 z(1)^2 - 1}, \quad \omega_{k+1} = \frac{1}{1 - \frac{\omega_k}{4 z(1)^2}} \quad (k \ge 2)$$
3. **Actualización vectorial de tres términos:**$$u^{(k+1)} = \omega_{k+1} \Big[ \gamma(G u^{(k)} + c) + (1 - \gamma)u^{(k)} - u^{(k-1)} \Big] + u^{(k-1)}$$Este algoritmo tiene exactamente el mismo costo de multiplicación matriz-vector por paso que el método base, requiere solo un vector adicional en memoria ($u^{(k-1)}$), y reduce el factor de convergencia asintótico de $\rho \approx 1 - \epsilon$ a $\approx 1 - \sqrt{2\epsilon}$.
# 2. Optimización y Métodos de Gradiente (Kincaid 4.7 & Strang VI.1, VI.4)
Cuando la matriz $A \in \mathbb{R}^{n \times n}$ es **simétrica y definida positiva (SPD)** ($A^T = A$ y $x^T Ax > 0$ para todo $x \neq 0$), abandonamos la perspectiva estática de punto fijo y reformulamos la resolución de $Ax = b$ como un problema variacional de optimización convexa sin restricciones.
## 2.1. Equivalencia Variacional: Sistema Lineal $\iff$ Minimización Cuadrática
Definimos el funcional cuadrático de energía $\phi: \mathbb{R}^n \to \mathbb{R}$:
$$\phi(x) = \frac{1}{2} x^T A x - b^T x$$
Calculando su vector gradiente y su matriz Hessiana (Kincaid 4.7, Strang I.7 y VI.1):
$$\nabla \phi(x) = \frac{1}{2}(A + A^T)x - b \underset{(A=A^T)}{=} Ax - b = -r(x), \qquad \nabla^2 \phi(x) = A \succ 0$$
>[!note] **Teorema de Equivalencia Variacional:**
> Si $A$ es SPD, el vector $x^* \in \mathbb{R}^n$ es solución de $Ax^* = b$ si y solo si $x^*$ es el único minimizador global estricto de $\phi(x)$.

_Demostración algebraica directa:_ Para cualquier $x \in \mathbb{R}^n$, sea $e = x - x^*$ el vector error. Sustituyendo $b = Ax^*$ en $\phi(x)$ y usando que $A^T = A$:
$$\begin{aligned} \phi(x) &= \frac{1}{2} x^T A x - (Ax^*)^T x = \frac{1}{2} x^T A x - (x^*)^T A x \\ &= \frac{1}{2}(x - x^*)^T A (x - x^*) - \frac{1}{2}(x^*)^T A x^* = \frac{1}{2} |x - x^*|_A^2 + \phi(x^*) \end{aligned}$$
De aquí se desprende la identidad fundamental:
$$\phi(x) - \phi(x^*) = \frac{1}{2} (x - x^*)^T A (x - x^*) = \frac{1}{2} |e|_A^2 \ge 0$$
Por ser $A$ definida positiva, $\|e\|_A^2 > 0$ para todo $x \neq x^*$, alcanzándose el mínimo absoluto $\phi(x^*) = -\frac{1}{2}(x^*)^T A x^*$ única y exclusivamente cuando $e = 0 \iff x = x^*$. Además, minimizar $\phi(x)$ en cualquier subespacio es exactamente equivalente a minimizar la distancia a la solución medida en la norma de energía $\|x - x^*\|_A$. $\blacksquare$
## 2.2. Método de Descenso por Gradiente (Steepest Descent — Máximo Descenso)
En la iteración $k$, situados en $x^{(k)}$, la dirección unitaria en la que $\phi$ decrece más rápidamente es la opuesta al gradiente:
$$-\nabla \phi(x^{(k)}) = b - Ax^{(k)} = r^{(k)}$$
Proponemos la actualización unidimensional $x^{(k+1)} = x^{(k)} + \alpha_k r^{(k)}$.
### Derivación Analítica del Tamaño de Paso Exacto (Exact Line Search)
Buscamos el escalar $\alpha_k \in \mathbb{R}$ que minimiza la función univariada $g(\alpha) = \phi(x^{(k)} + \alpha r^{(k)})$. Expandiendo cuadráticamente (con $r = r^{(k)}$ por brevedad):
$$\begin{aligned} g(\alpha) &= \frac{1}{2}(x^{(k)} + \alpha r)^T A (x^{(k)} + \alpha r) - b^T (x^{(k)} + \alpha r) \\ &= \left( \frac{1}{2} x^{(k)T} A x^{(k)} - b^T x^{(k)} \right) + \alpha r^T (Ax^{(k)} - b) + \frac{1}{2} \alpha^2 r^T A r \\ &= \phi(x^{(k)}) - \alpha r^T r + \frac{1}{2} \alpha^2 r^T A r \end{aligned}$$
Como $A$ es SPD y $r \neq 0$, el coeficiente cuadrático es estrictamente positivo ($r^T A r > 0$), por lo que $g(\alpha)$ es una parábola convexa hacia arriba. Imponiendo la condición de estacionariedad $g'(\alpha_k) = 0$:
$$\frac{d}{d\alpha} \phi(x^{(k)} + \alpha r^{(k)}) \Big|_{\alpha = \alpha_k} = -r^{(k)T} r^{(k)} + \alpha_k r^{(k)T} A r^{(k)} = 0 \implies \boxed{\alpha_k = \frac{r^{(k)T} r^{(k)}}{r^{(k)T} A r^{(k)}}}$$
Las fórmulas recursivas del método de Descenso por Gradiente son:
$$x^{(k+1)} = x^{(k)} + \alpha_k r^{(k)}, \qquad r^{(k+1)} = b - Ax^{(k+1)} = r^{(k)} - \alpha_k A r^{(k)}$$
### Ortogonalidad de Residuos Consecutivos y Geometría del Zigzag
Calculando el producto interno entre dos residuos sucesivos:
$$\langle r^{(k+1)}, r^{(k)} \rangle = r^{(k+1)T} r^{(k)} = (r^{(k)} - \alpha_k A r^{(k)})^T r^{(k)} = r^{(k)T} r^{(k)} - \left(\frac{r^{(k)T} r^{(k)}}{r^{(k)T} A r^{(k)}}\right) r^{(k)T} A r^{(k)} = 0$$
Geométricamente, las curvas de nivel $\{x : \phi(x) = c\}$ son hiperelipsoides concéntricos centrados en $x^*$cuyos ejes principales están alineados con los autovectores de $A$, y cuyas longitudes de semiejes son proporcionales a $1/\sqrt{\lambda_i}$ (Strang I.7 y VI.4). El paso óptimo avanza a lo largo de la normal a la curva de nivel en $x^{(k)}$ hasta el punto donde la recta de búsqueda es **tangente** a una nueva curva de nivel inferior. En ese punto de tangencia, el nuevo gradiente $\nabla \phi(x^{(k+1)}) = -r^{(k+1)}$ es forzosamente ortogonal a la dirección de búsqueda previa $r^{(k)}$.
Cuando el número de condición $\kappa(A) = \lambda_{\max}/\lambda_{\min} \gg 1$, los elipsoides son valles estrechos y muy alargados. A menos que $r^{(0)}$ apunte exactamente a lo largo de un eje principal (autovector), los pasos sucesivos rebotan perpendicularmente de una pared del valle a la otra, generando una **trayectoria en zigzag** patológicamente lenta (Strang VI.4).
### Cota de Convergencia y Desigualdad de Kantorovich
Usando la relación $\phi(x^{(k)}) - \phi(x^*) = \frac{1}{2}\|e^{(k)}\|_A^2$ y la expresión de $g(\alpha_k)$evaluada en el paso óptimo $\alpha_k$:
$$\frac{1}{2}\|e^{(k+1)}\|_A^2 - \frac{1}{2}\|e^{(k)}\|_A^2 = \phi(x^{(k+1)}) - \phi(x^{(k)}) = -\frac{1}{2} \frac{(r^{(k)T} r^{(k)})^2}{r^{(k)T} A r^{(k)}}$$
Observando que $\|e^{(k)}\|_A^2 = e^{(k)T} A e^{(k)} = (-A^{-1}r^{(k)})^T A (-A^{-1}r^{(k)}) = r^{(k)T} A^{-1} r^{(k)}$, dividimos por $\frac{1}{2}\|e^{(k)}\|_A^2$ para obtener la reducción exacta del error en un paso:
$$\|e^{(k+1)}\|_A^2 = \|e^{(k)}\|_A^2 \left[ 1 - \frac{(r^{(k)T} r^{(k)})^2}{(r^{(k)T} A r^{(k)})(r^{(k)T} A^{-1} r^{(k)})} \right]$$
Invocamos ahora la **Desigualdad de Kantorovich** (Kincaid 4.7): para toda matriz $A$ SPD con autovalores $0 < \lambda_{\min} = \lambda_1 \le \dots \le \lambda_n = \lambda_{\max}$ y todo vector $r \neq 0$:
$$\frac{(r^T A r)(r^T A^{-1} r)}{(r^T r)^2} \le \frac{(\lambda_{\max} + \lambda_{\min})^2}{4 \lambda_{\max} \lambda_{\min}} = \frac{(\kappa + 1)^2}{4\kappa}, \qquad \text{donde } \kappa = \kappa_2(A) = \frac{\lambda_{\max}}{\lambda_{\min}}$$
Invirtiendo la desigualdad y sustituyendo en la ecuación del error:
$$\|e^{(k+1)}\|_A^2 \le \|e^{(k)}\|_A^2 \left[ 1 - \frac{4\kappa}{(\kappa + 1)^2} \right] = \|e^{(k)}\|_A^2 \left( \frac{\kappa - 1}{\kappa + 1} \right)^2$$
Tomando raíz cuadrada e iterando desde $k=0$, obtenemos la cota clásica del Descenso por Gradiente:
$$\|e^{(k)}\|_A \le \left( \frac{\kappa - 1}{\kappa + 1} \right)^k \|e^{(0)}\|_A$$
_(Nota de conexión con Strang VI.4: Si en lugar de búsqueda lineal exacta se usa un paso fijo constante $s$, la matriz de iteración sobre el error es $I - sA$. El paso fijo óptimo que minimiza $\max_i |1 - s\lambda_i|$ es $s_{\text{opt}} = \frac{2}{\lambda_{\max} + \lambda_{\min}}$, que arroja exactamente el mismo factor de contracción en norma 2: $\frac{\lambda_{\max} - \lambda_{\min}}{\lambda_{\max} + \lambda_{\min}} = \frac{\kappa - 1}{\kappa + 1}$)._
## 2.3. Método de Gradiente Conjugado (CG — Hestenes & Stiefel, 1952)
Para eliminar de raíz el defecto del zigzag (donde optimizar en una nueva dirección arruina la optimalidad alcanzada en direcciones previas), el Método de Gradiente Conjugado exige que las direcciones de búsqueda $\{p_0, p_1, \dots, p_{k-1}\}$ sean mutuamente **$A$-ortogonales** (o **$A$-conjugadas**) (Kincaid 4.7, Strang II.1 y VI.4):
$$\langle p_i, p_j \rangle_A = p_i^T A p_j = 0 \qquad \forall i \neq j$$
### Subespacios de Krylov y Construcción del Algoritmo
Definimos el **$k$-ésimo Subespacio de Krylov** generado por $A$ y el residuo inicial $r^{(0)}$ como:
$$\mathcal{K}_k(A, r^{(0)}) = \text{span}\{r^{(0)}, A r^{(0)}, A^2 r^{(0)}, \dots, A^{k-1} r^{(0)}\}$$
En cada paso $k$, CG actualiza la solución en la dirección $p_k$: $x^{(k+1)} = x^{(k)} + \alpha_k p_k$, eligiendo $\alpha_k$ para minimizar $\phi(x^{(k)} + \alpha p_k)$. Derivando respecto de $\alpha$:
$$\frac{d}{d\alpha}\phi(x^{(k)} + \alpha p_k) = p_k^T (Ax^{(k)} - b) + \alpha p_k^T A p_k = -p_k^T r^{(k)} + \alpha p_k^T A p_k = 0 \implies \alpha_k = \frac{p_k^T r^{(k)}}{p_k^T A p_k}$$
La siguiente dirección $p_{k+1}$ se construye aplicando $A$-ortogonalización de Gram-Schmidt al nuevo residuo $r^{(k+1)}$ contra la dirección previa $p_k$:
$$p_{k+1} = r^{(k+1)} + \beta_k p_k$$
Imponiendo la condición de $A$-conjugación $p_{k+1}^T A p_k = 0$, despejamos $\beta_k = -\frac{r^{(k+1)T} A p_k}{p_k^T A p_k}$.
La magia algebraica de CG (derivada de que $A$ es simétrica y su proyección en el espacio de Krylov es tridiagonal, como en la iteración de Lanczos de Strang II.1) es que **conjugar únicamente contra la última dirección $p_k$ hace que $p_{k+1}$ resulte automáticamente $A$-conjugada a todas las direcciones anteriores $p_0, \dots, p_{k-1}$**.
>[!note] **Teorema de Propiedades Algebraicas y Terminación Finita de CG (Kincaid 4.7):**
> En aritmética exacta, mientras $r^{(k)} \neq 0$, las sucesiones generadas por CG satisfacen para todo $k \ge 1$:
> 1. **Igualdad de subespacios:** $\text{span}\{r^{(0)}, \dots, r^{(k-1)}\} = \text{span}\{p_0, \dots, p_{k-1}\} = \mathcal{K}_k(A, r^{(0)})$.
> 2. **Ortogonalidad de residuos con el subespacio previo:** $r^{(k)T} p_j = 0$ para todo $j < k$.
> 3. **Ortogonalidad mutua de todos los residuos:** $r^{(k)T} r^{(j)} = 0$ para todo $j < k$.
> 4. **$A$-conjugación mutua de todas las direcciones:** $p_k^T A p_j = 0$ para todo $j < k$.
> 5. **Fórmulas simplificadas de un solo producto matriz-vector:**
>$$\alpha_k = \frac{r^{(k)T} r^{(k)}}{p_k^T A p_k}, \qquad \beta_k = \frac{r^{(k+1)T} r^{(k+1)}}{r^{(k)T} r^{(k)}}$$
> 6. **Terminación finita:** El algoritmo converge a la solución exacta $x^*$ en a lo sumo $n$ iteraciones (más aún, en a lo sumo $m \le n$ pasos, donde $m$ es el número de autovalores distintos de $A$).
### Demostración Completa (Inducción y Terminación)
- **Paso base ($k=1$):** Con $p_0 = r^{(0)}$, tenemos $r^{(1)} = r^{(0)} - \alpha_0 A p_0 \in \mathcal{K}_2(A, r^{(0)})$. Por definición de $\alpha_0$, $p_0^T r^{(1)} = r^{(0)T} r^{(1)} = 0$. Además, por definición de $\beta_0$, $p_1^T A p_0 = 0$.
 - **Paso inductivo:** Supongamos ciertas las propiedades hasta el índice $k$.
 - _Prueba de (2) y (3) para $k+1$:_ Como $r^{(k+1)} = r^{(k)} - \alpha_k A p_k$, para $j = k$ tenemos $p_k^T r^{(k+1)} = p_k^T r^{(k)} - \alpha_k p_k^T A p_k = 0$ por la elección de $\alpha_k$. Para $j < k$, multiplicando por $p_j^T$:$$p_j^T r^{(k+1)} = \underbrace{p_j^T r^{(k)}}_{0 \text{ (por H.I.)}} - \alpha_k \underbrace{p_j^T A p_k}_{0 \text{ (por H.I.)}} = 0$$Como $\text{span}\{r^{(0)}, \dots, r^{(k)}\} = \text{span}\{p_0, \dots, p_k\}$, el hecho de que $r^{(k+1)} \perp p_j$ para todo $j \le k$ implica inmediatamente que $r^{(k+1)} \perp r^{(j)}$ para todo $j \le k$.
- _Prueba de (5):_ Como $p_k = r^{(k)} + \beta_{k-1} p_{k-1}$ y $r^{(k)T} p_{k-1} = 0$, resulta $p_k^T r^{(k)} = r^{(k)T} r^{(k)}$, lo que demuestra $\alpha_k = \frac{r^{(k)T} r^{(k)}}{p_k^T A p_k}$. A su vez, despejando $A p_k = \frac{1}{\alpha_k}(r^{(k)} - r^{(k+1)})$ y multiplicando por $r^{(k+1)T}$ (usando que $r^{(k+1)T} r^{(k)} = 0$):
$$r^{(k+1)T} A p_k = -\frac{1}{\alpha_k} r^{(k+1)T} r^{(k+1)} \implies \beta_k = -\frac{r^{(k+1)T} A p_k}{p_k^T A p_k} = \frac{r^{(k+1)T} r^{(k+1)}}{\alpha_k p_k^T A p_k} = \frac{r^{(k+1)T} r^{(k+1)}}{r^{(k)T} r^{(k)}}$$
- _Prueba de (4) para $k+1$:_ Para $j = k$, $p_{k+1}^T A p_k = 0$ por construcción de $\beta_k$. Para $j < k$:$$p_{k+1}^T A p_j = r^{(k+1)T} A p_j + \beta_k \underbrace{p_k^T A p_j}_{0 \text{ (por H.I.)}} = r^{(k+1)T} \left( \frac{r^{(j)} - r^{(j+1)}}{\alpha_j} \right) = 0$$ya que para $j < k$ se cumple $j+1 \le k$, y ya probamos que $r^{(k+1)}$ es ortogonal a todos los residuos de índice $\le k$.
- **Demostración de Terminación Finita:** Los vectores residuos $\{r^{(0)}, r^{(1)}, \dots, r^{(k)}\}$ forman un conjunto de vectores mutuamente ortogonales en $\mathbb{R}^n$. En un espacio vectorial de dimensión $n$ no pueden existir más de $n$ vectores no nulos mutuamente ortogonales (pues son linealmente independientes). Por lo tanto, necesariamente existe un índice $m \le n$ tal que $r^{(m)} = 0 \implies Ax^{(m)} = b \implies x^{(m)} = x^*$. $\blacksquare$
### Optimalidad en Espacios de Krylov y Cota Multiplicativa de Convergencia con $\sqrt{\kappa}$
Dado que $x^{(k)} \in x^{(0)} + \mathcal{K}_k(A, r^{(0)})$ y el error $e^{(k)} = x^{(k)} - x^*$ es $A$-ortogonal a todo el espacio $\mathcal{K}_k(A, r^{(0)})$ (pues $\langle e^{(k)}, p_j \rangle_A = -r^{(k)T} p_j = 0$ para $j < k$), por el teorema de proyección ortogonal en el espacio de Hilbert $(\mathbb{R}^n, \langle \cdot, \cdot \rangle_A)$, $x^{(k)}$ es la mejor aproximación a $x^*$ en norma de energía:
$$\|e^{(k)}\|_A = \min_{x \in x^{(0)} + \mathcal{K}_k(A, r^{(0)})} \|x - x^*\|_A = \min_{\substack{P_k \in \Pi_k \\ P_k(0) = 1}} \|P_k(A) e^{(0)}\|_A$$
donde hemos usado que todo vector en $x^{(0)} + \mathcal{K}_k(A, r^{(0)})$ tiene un error de la forma $e^{(0)} + c_1 A e^{(0)} + \dots + c_k A^k e^{(0)} = P_k(A) e^{(0)}$ con $P_k(0) = 1$ (ya que $r^{(0)} = -Ae^{(0)}$).
Expandiendo $e^{(0)}$ en la base ortonormal de autovectores de $A$, obtenemos la cota espectral:
$$\|e^{(k)}\|_A \le \left( \min_{\substack{P_k \in \Pi_k \\ P_k(0) = 1}} \max_{\lambda \in [\lambda_{\min}, \lambda_{\max}]} |P_k(\lambda)| \right) \|e^{(0)}\|_A$$
(Observación pedagógica: Si $A$ tiene solo $m < n$ autovalores distintos $\lambda_1, \dots, \lambda_m$, el polinomio interpolador de Lagrange $P_m(t) = \prod_{i=1}^m (1 - t/\lambda_i)$ tiene grado $m$, cumple $P_m(0)=1$ y se anula en todos los autovalores de $A$, probando que $\|e^{(m)}\|_A = 0$ en a lo sumo $m$ pasos).
Para acotar el problema minimax en todo el intervalo continuo $[\lambda_{\min}, \lambda_{\max}]$,empleamos nuevamente el polinomio de Chebyshev trasladado $P_k(t) = \frac{T_k\left(\frac{\lambda_{\max} + \lambda_{\min} - 2t}{\lambda_{\max} - \lambda_{\min}}\right)}{T_k\left(\frac{\lambda_{\max} + \lambda_{\min}}{\lambda_{\max} - \lambda_{\min}}\right)}$. Evaluando el denominador en $z_0 = \frac{\lambda_{\max} + \lambda_{\min}}{\lambda_{\max} - \lambda_{\min}} = \frac{\kappa + 1}{\kappa - 1} > 1$ mediante la identidad explícita $T_k(z) = \frac{1}{2}\big[ (z + \sqrt{z^2 - 1})^k + (z - \sqrt{z^2 - 1})^k \big]$:
$$z_0 + \sqrt{z_0^2 - 1} = \frac{\kappa + 1 + \sqrt{(\kappa+1)^2 - (\kappa-1)^2}}{\kappa - 1} = \frac{\kappa + 1 + 2\sqrt{\kappa}}{\kappa - 1} = \frac{(\sqrt{\kappa} + 1)^2}{(\sqrt{\kappa} - 1)(\sqrt{\kappa} + 1)} = \frac{\sqrt{\kappa} + 1}{\sqrt{\kappa} - 1}$$
Por consiguiente, $T_k(z_0) > \frac{1}{2} \left( \frac{\sqrt{\kappa} + 1}{\sqrt{\kappa} - 1} \right)^k$, lo que demuestra la célebre **cota de convergencia del Gradiente Conjugado** (Kincaid 4.7, Strang II.1 ec. 8):
$$\boxed{\|e^{(k)}\|_A \le 2 \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^k \|e^{(0)}\|_A, \qquad \kappa = \kappa_2(A) = \frac{\lambda_{\max}}{\lambda_{\min}}}$$
**Impacto cuantitativo ($\kappa$ vs. $\sqrt{\kappa}$):** Si $\kappa(A) = 10^4$, el factor de contracción de Descenso por Gradiente es $\frac{9999}{10001} \approx 0.9998$ (requiere $\sim 23\,000$ iteraciones para ganar 2 dígitos), mientras que en CG el factor es $\frac{99}{101} \approx 0.9802$ (requiere apenas $\sim 230$ iteraciones, 100 veces más rápido con el mismo costo por paso).
## 2.4. Precondicionamiento (PCG — Preconditioned Conjugate Gradient)
Cuando $\kappa(A)$ es extremadamente grande, incluso $\sqrt{\kappa(A)}$ puede degradar la velocidad de CG. Como explica Strang (Sección II.1) y Kincaid (Sección 4.7), la estrategia de **precondicionamiento** introduce una matriz simétrica y definida positiva $M$ (llamada _precondicionador_, notada $P$ en Strang II.1) que cumple dos requisitos simultáneos:
1. Los sistemas lineales $Mz = r$ son muy rápidos de resolver ($\mathcal{O}(n)$ flops).
2. $M$ aproxima espectralmente a $A$, de modo que $\kappa(M^{-1}A) \ll \kappa(A)$ (o los autovalores de $M^{-1}A$ quedan fuertemente agrupados alrededor de 1).
Para preservar la simetría indispensable en CG, no aplicamos CG directamente a $M^{-1}Ax = M^{-1}b$ (pues $M^{-1}A$ no suele ser simétrica en el producto interno euclídeo, como advierte Strang en I.10 y II.1), sino al sistema simétrico transformado por congruencia con $M^{-1/2}$ (o factor de Cholesky $M = CC^T$):
$$\underbrace{(M^{-1/2} A M^{-1/2})}_{\tilde{A}} \underbrace{(M^{1/2} x)}_{\tilde{x}} = \underbrace{M^{-1/2} b}_{\tilde{b}}$$
Al aplicar CG sobre $\tilde{A}\tilde{x} = \tilde{b}$ y regresar a las variables originales ($z^{(k)} = M^{-1}r^{(k)}$), los factores $M^{-1/2}$ se cancelan algebraicamente. **Nunca necesitamos calcular $M^{1/2}$**: en cada iteración solo resolvemos un sistema auxiliar $Mz^{(k)} = r^{(k)}$.
### Algoritmo PCG (En variables originales)
- Inicializar $x^{(0)}$, $r^{(0)} = b - Ax^{(0)}$, resolver $Mz^{(0)} = r^{(0)}$, hacer $p_0 = z^{(0)}$.
 - Para $k = 0, 1, \dots$:
 $$\alpha_k = \frac{r^{(k)T} z^{(k)}}{p_k^T A p_k}, \quad x^{(k+1)} = x^{(k)} + \alpha_k p_k, \quad r^{(k+1)} = r^{(k)} - \alpha_k A p_k$$
 $$\text{Resolver } M z^{(k+1)} = r^{(k+1)}, \quad \beta_k = \frac{r^{(k+1)T} z^{(k+1)}}{r^{(k)T} z^{(k)}}, \quad p_{k+1} = z^{(k+1)} + \beta_k p_k$$
### Estrategias Principales de Precondicionadores (Strang II.1 & Kincaid 4.7)
1. **Precondicionador de Jacobi (Diagonal):** $M = \text{diag}(a_{11}, \dots, a_{nn})$. Costo casi nulo ($\mathcal{O}(n)$ divisiones), totalmente paralelizable. Elimina el mal condicionamiento causado por escalas dispares entre variables (equilibrado diagonal).
2. **Precondicionador SSOR:** $M = Q_{\text{SSOR}} = \frac{1}{\omega(2-\omega)}(D + \omega L)D^{-1}(D + \omega L^T)$. No requiere ensamblaje adicional; resolver $Mz = r$ equivale a una sustitución hacia adelante y otra hacia atrás usando los propios elementos de $A$.
 3. **Factorización Incompleta LU / Cholesky ($\text{ILU}(0)$ / $\text{IC}(0)$):** Se ejecuta la eliminación gaussiana o de Cholesky sobre $A$, pero se fuerza a cero cualquier elemento en posiciones $(i,j)$ donde la matriz original $A$ tenía un cero ($a_{ij} = 0$). Así, $M = L_0 L_0^T$ (o $L_0 U_0$) preserva exactamente el mismo patrón de raleza (_sparsity pattern_) que $A$, evitando el _fill-in_ en memoria y aproximando muy bien la inversa de $A$.
 4. **Precondicionador Multimalla (Multigrid — Kincaid 9.8 & Strang II.1):** Los métodos iterativos clásicos (como Jacobi ponderado o Gauss-Seidel) atenúan rápidamente las componentes de error de **alta frecuencia** (oscilatorias en la malla), pero son muy lentos en las componentes de **baja frecuencia**(suaves). El método Multimalla proyecta (restringe) el residuo suavizado a una secuencia de mallas progresivamente más gruesas (donde los modos suaves pasan a ser oscilatorios y los sistemas son diminutos), resuelve allí la ecuación del defecto y prolonga (interpola) la corrección de vuelta a la malla fina. Como precondicionador, logra un número de condición acotado $\kappa(M^{-1}A) = \mathcal{O}(1)$independiente del tamaño de malla $h$, alcanzando complejidad óptima $\mathcal{O}(n)$.
# 3. Variantes Modernas y Aplicaciones en Ciencia de Datos y ML (Strang VI.4, VI.5, III.4)

En Ciencia de Datos y Aprendizaje Profundo (_Deep Learning_), buscamos minimizar una **función de pérdida**(o riesgo empírico) $L(x) = \frac{1}{N} \sum_{i=1}^N \ell_i(x)$ donde $x \in \mathbb{R}^d$ contiene millones de parámetros o pesos de una red neuronal y $N$ es la cantidad de muestras de entrenamiento (Strang, Prefacio y Parte VI). A diferencia del problema lineal clásico, aquí $L(x)$ no es cuadrática (a menudo ni siquiera es convexa) y calcular o invertir la matriz Hessiana $\nabla^2 L(x)$ de dimensión $d \times d$ es computacionalmente inviable. Nos apoyamos exclusivamente en métodos de primer orden basados en el gradiente $\nabla L(x)$.
## 3.1. Incorporación de Momento: Bola Pesada de Polyak y Aceleración de Nesterov

Para atenuar el efecto zigzag del descenso por gradiente en valles estrechos ($\kappa \gg 1$) sin calcular pasos exactos ni productos matriz-vector con la Hessiana, se introduce un término físico de **inercia o momento** (_momentum_) que acumula una media móvil exponencial de los gradientes pasados (Strang,Sección VI.4).
### 1. Método de la Bola Pesada de Polyak (Heavy Ball / Classical Momentum)
Imaginemos una esfera con masa rodando por la superficie del gráfico de $L(x)$:
$$x^{(k+1)} = x^{(k)} - s \nabla L(x^{(k)}) + \beta (x^{(k)} - x^{(k-1)})$$
Equivalente al sistema de dos pasos con un vector de velocidad (o término de actualización) $z^{(k)}$:
$$z^{(k+1)} = \beta z^{(k)} - s \nabla L(x^{(k)}), \qquad x^{(k+1)} = x^{(k)} + z^{(k+1)}$$
- **Mecanismo geométrico contra el zigzag:** En las direcciones transversales al valle (autovalores grandes de la Hessiana, $\lambda_{\max}$), el gradiente cambia de signo en cada iteración; al sumarse en $z^{(k)}$, las componentes opuestas se **cancelan mutuamente**, amortiguando las oscilaciones. En cambio, a lo largo del fondo del valle (autovalores pequeños, $\lambda_{\min}$), el gradiente apunta consistentemente hacia el mínimo; el término $\beta z^{(k)}$ **refuerza y acelera** el avance en esa dirección.
 - **Análisis espectral en el caso cuadrático ($L(x) = \frac{1}{2}x^T A x - b^T x$):** Como la actualización depende de $x^{(k)}$ y $x^{(k-1)}$, es una ecuación en diferencias de segundo orden (Kincaid 1.3, Strang VI.4). El estado ampliado de errores satisface:$$\begin{pmatrix} e^{(k+1)} \\ e^{(k)} \end{pmatrix} = \underbrace{\begin{pmatrix} (1+\beta)I - sA & -\beta I \\ I & 0 \end{pmatrix}}_{R_{\text{HB}}} \begin{pmatrix} e^{(k)} \\ e^{(k-1)} \end{pmatrix}$$Para cada autovalor $\lambda_i$ de $A$, los autovalores $r$ del bloque $2 \times 2$ correspondiente satisfacen la ecuación característica $r^2 - (1 + \beta - s\lambda_i)r + \beta = 0$. Eligiendo los parámetros óptimos que igualan el discriminante a cero en los extremos del espectro (haciendo que todas las raíces sean complejas conjugadas de módulo exactamente $\sqrt{\beta}$):$$s_{\text{opt}} = \left( \frac{2}{\sqrt{\lambda_{\max}} + \sqrt{\lambda_{\min}}} \right)^2, \qquad \beta_{\text{opt}} = \left( \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1} \right)^2 \implies \rho(R_{\text{HB}}) = \sqrt{\beta_{\text{opt}}} = \frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}$$¡La Bola Pesada alcanza asintóticamente la misma tasa óptima $\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}$que el Gradiente Conjugado y Chebyshev, utilizando coeficientes constantes!
### 2. Gradiente Acelerado de Nesterov (NAG — Nesterov Accelerated Gradient)
En el método de la Bola Pesada, la evaluación del gradiente ocurre en la posición actual $x^{(k)}$ _antes_ de aplicar el impulso de momento. Nesterov (1983) propuso una modificación sutil pero profunda (Strang VI.4):primero dar el salto inercial hasta un punto de **mirada adelantada** (_look-ahead point_) $y^{(k)}$, y evaluar el gradiente allí:
$$\begin{aligned} y^{(k)} &= x^{(k)} + \beta_k (x^{(k)} - x^{(k-1)}) \\ x^{(k+1)} &= y^{(k)} - s \nabla L(y^{(k)}) \end{aligned}$$
- **Ventaja física y teórica:** Si el impulso $\beta_k (x^{(k)} - x^{(k-1)})$ está por llevarnos cuesta arriba al otro lado del valle, Nesterov mide el gradiente en $y^{(k)}$ (donde la pendiente ya apunta de regreso hacia abajo) y frena o corrige el rumbo _antes_ de completar el paso. Mientras que la Bola Pesada solo garantiza su tasa óptima para funciones cuadráticas (y puede divergir en funciones convexas generales),el método de Nesterov logra la tasa óptima global de convergencia para funciones convexas suaves ($\mathcal{O}(1/k^2)$ en lugar de $\mathcal{O}(1/k)$) y $\left(1 - \frac{1}{\sqrt{\kappa}}\right)^k$ para funciones fuertemente convexas.
## 3.2. Descenso por Gradiente Estocástico (SGD) y Algoritmo de Kaczmarz
Cuando el conjunto de datos consta de $N \sim 10^6$ o $10^9$ muestras, la función de costo total es una suma finita masiva:
$$L(x) = \frac{1}{N} \sum_{i=1}^N \ell_i(x)$$
En el Descenso por Gradiente de lote completo (_Batch GD_), calcular un solo vector gradiente $\nabla L(x) = \frac{1}{N}\sum_{i=1}^N \nabla \ell_i(x)$ exige evaluar la red sobre los $N$ datos antes de dar un único paso. Esto es computacionalmente prohibitivo y redundante.
### 1. Formulación de SGD y Uso de Minilotes (Mini-batch SGD — Strang VI.5)
En cada iteración $k$, seleccionamos al azar (uniformemente o sin reemplazo en cada _época_) una única muestra $i_k \in \{1, \dots, N\}$ o un **minilote** (_mini-batch_) $B_k \subset \{1, \dots, N\}$ de tamaño moderado $b = |B_k|$ (típicamente $b \in \{32, 64, 128, 256\}$), y actualizamos los pesos usando el estimador estocástico insesgado del gradiente:
$$x^{(k+1)} = x^{(k)} - s_k g^{(k)}, \qquad \text{donde } g^{(k)} = \frac{1}{b} \sum_{i \in B_k} \nabla \ell_i(x^{(k)})$$
Como $\mathbb{E}[g^{(k)} \mid x^{(k)}] = \nabla L(x^{(k)})$, en promedio el paso desciende en la dirección correcta, pero con un costo computacional reducido en un factor $N/b$. Además, evaluar un minilote de $b$vectores en simultáneo convierte los productos matriz-vector $Ax$ en productos matriz-matriz $AX_B$,aprovechando al máximo el paralelismo de las unidades tensoriales en GPUs.
### 2. Análisis de Semiconvergencia (Fase de Transición vs. Fase de Fluctuación)
En el mínimo global $x^*$ de la pérdida total se cumple $\nabla L(x^*) = \frac{1}{N}\sum_{i=1}^N \nabla \ell_i(x^*) = 0$, pero **los gradientes individuales no se anulan**: $\nabla \ell_i(x^*) \neq 0$ (ya que ninguna curva ajusta exactamente todos los datos ruidosos a la vez). Sea $\sigma^2 = \mathbb{E}[\|\nabla \ell_i(x^*)\|^2]$ la varianza intrínseca del gradiente en el óptimo. Si mantenemos una tasa de aprendizaje (_learning rate_) constante $s_k = s$, el error cuadrático medio esperado para una función $\mu$-fuertemente convexa cumple una desigualdad del tipo (Strang VI.5):
$$\mathbb{E}\big[\|x^{(k)} - x^*\|^2\big] \le (1 - 2\mu s)^k \|x^{(0)} - x^*\|^2 + \frac{s \sigma^2}{2 \mu b}$$
Este resultado revela el comportamiento de **semiconvergencia** en dos fases:
1. **Fase 1 (Búsqueda rápida / Descenso transitorio):** Lejos de $x^*$, el término determinista $(1 - 2\mu s)^k \|e^{(0)}\|^2$ domina y el error cae exponencialmente rápido hacia la vecindad de la solución.
2. **Fase 2 (Piso de ruido / Bola de confusión):** Cerca de $x^*$, la varianza del muestreo domina y las iteraciones no convergen a un punto, sino que oscilan perpetuamente dentro de una "bola de ruido" de radio proporcional a $\sqrt{s/b}$. Para lograr convergencia exacta ($x^{(k)} \to x^*$), es necesario disminuir gradualmente el paso $s_k \to 0$ (condiciones de Robbins-Monro: $\sum_{k=1}^\infty s_k = \infty$ y $\sum_{k=1}^\infty s_k^2 < \infty$, por ejemplo $s_k = \frac{c}{k + k_0}$) o incrementar el tamaño de minilote $b$.
### 3. Equivalencia con el Algoritmo de Kaczmarz Aleatorizado (Strang II.1 y VI.5)
Apliquemos SGD con tamaño de lote $b=1$ al sistema lineal $Ax = b$, donde la $i$-ésima fila de $A$ es el vector $a_i^T$. Definimos la pérdida individual de la $i$-ésima ecuación normalizada por la norma de su fila:
$$\ell_i(x) = \frac{1}{2} \frac{(a_i^T x - b_i)^2}{\|a_i\|_2^2} \implies \nabla \ell_i(x) = \frac{a_i^T x - b_i}{\|a_i\|_2^2} a_i$$
Tomando un paso de SGD con tasa de aprendizaje unitaria $s_k = 1$ sobre la ecuación elegida $i = i_k$, obtenemos la célebre **Iteración de Kaczmarz** (Strang II.1 ec. 10 y VI.5):
$$x^{(k+1)} = x^{(k)} - \nabla \ell_i(x^{(k)}) = x^{(k)} + \frac{b_i - a_i^T x^{(k)}}{\|a_i\|_2^2} a_i$$
- **Interpretación geométrica (Proyecciones sucesivas en hiperplanos):** Cada ecuación $a_i^T x = b_i$define un hiperplano afín $H_i \subset \mathbb{R}^n$ cuyo vector normal es $a_i$. Se verifica inmediatamente que $a_i^T x^{(k+1)} = a_i^T x^{(k)} + (b_i - a_i^T x^{(k)}) = b_i$, y que el desplazamiento $x^{(k+1)} - x^{(k)}$ es paralelo a la normal $a_i$. Por lo tanto, **cada paso de Kaczmarz es la proyección ortogonal exacta de $x^{(k)}$ sobre el hiperplano $H_i$**. 
- **Teorema de Strohmer y Vershynin (2009 — Strang II.1):** Si en lugar de recorrer las filas cíclicamente ($i = 1, 2, \dots, m$), en cada paso se elige la fila $i$ al azar con **probabilidad proporcional a la norma al cuadrado de la fila**, $p_i = \frac{|a_i|_2^2}{|A|_F^2}$, el algoritmo de Kaczmarz Aleatorizado converge exponencialmente rápido en esperanza para cualquier sistema consistente:
$$\mathbb{E}\big[\|x^{(k)} - x^*\|_2^2\big] \le \left( 1 - \frac{1}{\kappa_{\text{dem}}(A)^2} \right)^k \|x^{(0)} - x^*\|_2^2, \qquad \text{donde } \kappa_{\text{dem}}(A) = \frac{\|A\|_F}{\sigma_{\min}(A)}$$
## 3.3. Métodos Adaptativos y Promediado de Pesos (Strang VI.5)

En problemas de alta dimensión con datos dispersos (como procesamiento de lenguaje natural o redes profundas), distintas componentes $x_j$ del vector de pesos tienen sensibilidades muy diferentes: algunas coordenadas reciben gradientes frecuentes y enormes, mientras que otras corresponden a características raras y reciben gradientes casi nulos. Usar un único escalar de paso $s_k$ para todas las coordenadas es subóptimo. Los **métodos adaptativos** escalan individualmente cada coordenada dividiendo por una estimación histórica de la magnitud de sus gradientes (actuando como un precondicionador diagonal dinámico).
### 1. ADAGRAD (Adaptive Gradient — Duchi, Hazan & Singer, 2011)
Sea $g^{(k)} = \nabla \ell_{B_k}(x^{(k)}) \in \mathbb{R}^d$ el gradiente estocástico en el paso $k$. ADAGRAD acumula la suma histórica de los cuadrados de los gradientes componente a componente (operación de Hadamard $\odot$):
$$G^{(k)} = G^{(k-1)} + g^{(k)} \odot g^{(k)} \quad \text{(es decir, } G_j^{(k)} = \sum_{\tau=1}^k \big(g_j^{(\tau)}\big)^2 \text{ para cada coordenada } j=1,\dots,d\text{)}$$
La regla de actualización modifica la tasa de aprendizaje de cada parámetro $j$ de forma inversamente proporcional a $\sqrt{G_j^{(k)}}$:
$$x_j^{(k+1)} = x_j^{(k)} - \frac{s}{\sqrt{G_j^{(k)}} + \epsilon} \, g_j^{(k)} \iff x^{(k+1)} = x^{(k)} - s \big( \text{diag}(G^{(k)})^{1/2} + \epsilon I \big)^{-1} g^{(k)}$$
donde $\epsilon \approx 10^{-8}$ evita la división por cero.
- **Ventaja y limitación:** Asigna pasos grandes a parámetros poco frecuentes y pasos pequeños a parámetros frecuentes. Sin embargo, como $G_j^{(k)}$ es una suma monótonamente creciente de términos positivos, en entrenamientos largos $G_j^{(k)} \to \infty$, haciendo que la tasa de aprendizaje efectiva colapse prematuramente a cero y detenga el aprendizaje.
### 2. ADAM (Adaptive Moment Estimation — Kingma & Ba, 2015)
Como destaca Strang en la Sección VI.5, **ADAM** es el optimizador estándar por excelencia en el entrenamiento de redes neuronales profundas. Resuelve el estancamiento de ADAGRAD reemplazando la suma acumulada por una **media móvil exponencial** de los gradientes al cuadrado (como en RMSProp), y la combina con la **media móvil exponencial de los gradientes** (momento de primer orden de Polyak):
1. **Estimación del primer momento (Media / Momentum, con $\beta_1 \approx 0.9$):**
$$m^{(k)} = \beta_1 m^{(k-1)} + (1 - \beta_1) g^{(k)}$$
2. **Estimación del segundo momento no centrado (Varianza / Escala, con $\beta_2 \approx 0.999$):**
$$v^{(k)} = \beta_2 v^{(k-1)} + (1 - \beta_2) \big(g^{(k)} \odot g^{(k)}\big)$$
3. **Corrección de sesgo inicial (Bias Correction):** Como se inicializa con $m^{(0)} = 0$ y $v^{(0)} = 0$, en las primeras iteraciones ambos vectores están fuertemente sesgados hacia el vector nulo (especialmente cuando $\beta_1, \beta_2$ son cercanos a 1). Dado que $\mathbb{E}[m^{(k)}] = (1 - \beta_1^k)\mathbb{E}[g]$, se corrige la escala dividiendo por $(1 - \beta^k)$:
$$\hat{m}^{(k)} = \frac{m^{(k)}}{1 - \beta_1^k}, \qquad \hat{v}^{(k)} = \frac{v^{(k)}}{1 - \beta_2^k}$$
 4. **Actualización adaptativa de pesos:**
 $$\boxed{x^{(k+1)} = x^{(k)} - \frac{s}{\sqrt{\hat{v}^{(k)}} + \epsilon} \odot \hat{m}^{(k)}}$$
### 3. Promediado de Iteraciones (Polyak-Ruppert Averaging y SWA)
Para eliminar las oscilaciones de la Fase 2 del descenso estocástico sin sacrificar la velocidad de exploración de un paso $s_k$ moderadamente grande, Polyak y Ruppert demostraron que el **promedio ergódico de los pesos visitados** a lo largo de la trayectoria:
$$\bar{x}^{(K)} = \frac{1}{K} \sum_{k=1}^K x^{(k)} \quad \text{o su actualización recursiva } \bar{x}^{(k+1)} = \frac{k}{k+1}\bar{x}^{(k)} + \frac{1}{k+1}x^{(k+1)}$$
reduce la varianza del ruido estocástico a velocidad óptima $\mathcal{O}(1/K)$. En Aprendizaje Profundo moderno, esta idea da lugar a **SWA (Stochastic Weight Averaging)** (Strang VI.5): durante el tramo final del entrenamiento, al promediar los pesos $x^{(k)}$ muestreados alrededor de un valle, SWA guía la solución hacia el centro de regiones más **planas y anchas** de la superficie de pérdida, lo cual mejora drásticamente la capacidad de generalización sobre datos de prueba no vistos.
## 3.4. Mínimos Cuadrados No Lineales y Módulos de Partición para Regularización $\ell_1$
### 1. Mínimos Cuadrados No Lineales y Algoritmo de Levenberg-Marquardt (Strang VI.1 & Kincaid 3.2)
Consideremos el ajuste de un modelo no lineal a datos mediante la minimización de la suma de residuos al cuadrado:
$$\min_{x \in \mathbb{R}^n} \phi(x) = \frac{1}{2} \|r(x)\|_2^2 = \frac{1}{2} \sum_{i=1}^m r_i(x)^2$$
donde $r: \mathbb{R}^n \to \mathbb{R}^m$ es el vector de residuos no lineales. Sea $J(x) \in \mathbb{R}^{m \times n}$ la **matriz Jacobiana** de $r(x)$, cuyas entradas son $J_{ij}(x) = \frac{\partial r_i}{\partial x_j}$. El gradiente y la matriz Hessiana exacta de $\phi(x)$ son:
$$\nabla \phi(x) = J(x)^T r(x), \qquad \nabla^2 \phi(x) = J(x)^T J(x) + \sum_{i=1}^m r_i(x) \nabla^2 r_i(x)$$
- **Método de Gauss-Newton:** Linealiza el residuo alrededor de $x^{(k)}$ mediante Taylor: $r(x^{(k)} + p) \approx r(x^{(k)}) + J_k p$ (con $J_k = J(x^{(k)})$), lo cual equivale a despreciar el término de segundas derivadas $\sum r_i \nabla^2 r_i$ en la Hessiana. El paso $p_{\text{GN}}$ resuelve las ecuaciones normales linealizadas:
$$(J_k^T J_k) p_{\text{GN}} = -J_k^T r(x^{(k)})$$
Si $J_k$ es de rango incompleto o mal condicionada, o si estamos lejos de la solución (residuos $r_i$grandes), $J_k^T J_k$ es casi singular y el paso de Gauss-Newton explota o diverge.
 - **Algoritmo de Levenberg-Marquardt (Amortiguamiento / Región de Confianza):** Introduce un parámetro de regularización (amortiguamiento de Tikhonov / _ridge_) $\lambda_k > 0$ sumado a la diagonal de $J_k^T J_k$ (Strang II.2 ec. 19 y VI.1):$$\big( J_k^T J_k + \lambda_k I \big) p_{\text{LM}} = -J_k^T r(x^{(k)})$$_(Variante de Marquardt: usar $(J_k^T J_k + \lambda_k \text{diag}(J_k^T J_k)) p_{\text{LM}} = -J_k^T r(x^{(k)})$ para invariancia de escala)._
- **Propiedad de interpolación continua:** La matriz $J_k^T J_k + \lambda_k I$ es siempre **simétrica y definida positiva** para todo $\lambda_k > 0$.
- Cuando $\lambda_k \to 0$ (cerca del mínimo, donde el modelo cuadrático es fiel), $p_{\text{LM}} \to p_{\text{GN}}$ y el método recupera la **convergencia rápida (casi cuadrática) de Gauss-Newton**.
- Cuando $\lambda_k \to \infty$ (lejos del mínimo o en regiones no convexas donde el paso falla),$p_{\text{LM}} \approx -\frac{1}{\lambda_k} J_k^T r(x^{(k)}) = -\frac{1}{\lambda_k} \nabla \phi(x^{(k)})$, comportándose como un paso corto y seguro de **Descenso por Gradiente**, lo que garantiza convergencia global.
### 2. Regularización $\ell_1$ (LASSO), Umbral Suave (Soft-Thresholding) y ADMM (Strang I.11, I.12, III.4)
En problemas inversos, _Compressed Sensing_ y selección de características estadísticas, buscamos soluciones **ralas** (_sparse_, con muchas componentes exactamente iguales a cero) resolviendo el problema **LASSO** (Tibshirani, 1996 — Strang I.12 ec. 6 y III.4):
$$\min_{x \in \mathbb{R}^n} \frac{1}{2} \|Ax - b\|_2^2 + \lambda \|x\|_1$$
Como explicó Strang en la Sección I.11 (Figura I.16), las bolas de nivel de la norma $\ell_1$ son politopos (diamantes con vértices puntiagudos sobre los ejes coordenados); al tocar el subespacio afín de los datos, el contacto ocurre en los vértices o aristas, anulando exactamente las coordenadas irrelevantes. Sin embargo, $\|x\|_1$ **no es diferenciable en cero**, por lo que el descenso por gradiente ordinario no puede aplicarse directamente.
- **Operador Proximal de $\ell_1$: Umbral Suave (Soft-Thresholding / Shrinkage — Strang III.4):** Consideremos primero el problema desacoplado elemental (cuando $A = I$): $\min_x \frac{1}{2}\|x - v\|_2^2 + \alpha \|x\|_1 = \sum_{i=1}^n \left( \frac{1}{2}(x_i - v_i)^2 + \alpha |x_i| \right)$. Usando subcliferenciales en cada coordenada $x_i$ ($0 \in x_i - v_i + \alpha \partial |x_i|$, donde $\partial |x_i| = \text{sign}(x_i)$ si $x_i \neq 0$ y $[-1, 1]$ si $x_i = 0$), la solución analítica cerrada está dada por el operador de **Umbral Suave** (_Soft-Thresholding_) $S_\alpha(v)$: $$\boxed{[S_\alpha(v)]_i = \text{sign}(v_i) \max(|v_i| - \alpha, 0) = \begin{cases} v_i - \alpha & \text{si } v_i > \alpha \\ 0 & \text{si } |v_i| \le \alpha \\ v_i + \alpha & \text{si } v_i < -\alpha \end{cases}}$$Este operador crea una "zona muerta" $[-\alpha, \alpha]$ que envía a cero exacto todas las componentes pequeñas (generando raleza) y contrae hacia el origen en una magnitud $\alpha$ las componentes que superan el umbral. Aplicar un paso de descenso por gradiente sobre la parte suave $\frac{1}{2}\|Ax-b\|_2^2$ seguido del operador $S_{s\lambda}$ da origen al algoritmo **ISTA** (_Iterative Shrinkage-Thresholding Algorithm_, o su versión acelerada con momento de Nesterov, **FISTA**):
$$x^{(k+1)} = S_{s\lambda} \Big( x^{(k)} - s A^T(Ax^{(k)} - b) \Big)$$ 
- **Método de Multiplicadores de Direcciones Alternadas (ADMM — Strang I.12 y III.4):** Para separar completamente la dificultad de la matriz $A$ (término cuadrático $\ell_2$) de la no diferenciabilidad de la norma $\ell_1$, el algoritmo **ADMM** aplica una **partición de variables** (_variable splitting_) reescribiendo el problema LASSO con dos vectores acoplados por una restricción lineal:$$\min_{x, z \in \mathbb{R}^n} \frac{1}{2}\|Ax - b\|_2^2 + \lambda \|z\|_1 \qquad \text{sujeto a } x - z = 0$$Formamos el **Lagrangiano Aumentado** con parámetro de penalización $\rho > 0$ y variable dual escalada $u = (1/\rho)y$: $$\mathcal{L}_\rho(x, z, u) = \frac{1}{2}\|Ax - b\|_2^2 + \lambda \|z\|_1 + \frac{\rho}{2}\|x - z + u\|_2^2 - \frac{\rho}{2}\|u\|_2^2$$ADMM minimiza $\mathcal{L}_\rho$ de manera alternada (un bloque a la vez, al estilo Gauss-Seidel) seguido de una actualización del multiplicador dual $u$:
	1. **Paso de minimización en $x$ (Sistema lineal SPD suave $\ell_2$):** $$x^{(k+1)} = \arg\min_x \left( \frac{1}{2}\|Ax - b\|_2^2 + \frac{\rho}{2}\|x - z^{(k)} + u^{(k)}\|_2^2 \right) = \big(A^T A + \rho I\big)^{-1} \Big( A^T b + \rho(z^{(k)} - u^{(k)}) \Big)$$_(Nota computacional: La matriz $A^T A + \rho I$ es SPD y constante en todas las iteraciones; se puede factorizar una sola vez por Cholesky al inicio o resolver con CG)._
	 2. **Paso de minimización en $z$ (Evaluación directa por Umbral Suave $\ell_1$):**
		$$z^{(k+1)} = \arg\min_z \left( \lambda \|z\|_1 + \frac{\rho}{2}\|x^{(k+1)} - z + u^{(k)}\|_2^2 \right) = S_{\lambda/\rho}\big(x^{(k+1)} + u^{(k)}\big)$$
	 3. **Actualización de la variable dual (Suma acumulada de residuos de acoplamiento):**
	    $$u^{(k+1)} = u^{(k)} + \big(x^{(k+1)} - z^{(k+1)}\big)$$
# 4. Síntesis Comparativa y Aplicaciones Prácticas
## 4.1. Tabla Comparativa Global de Métodos
_(Convenciones: $n$ = dimensión del vector de incógnitas $x$; $\text{nnz}(A)$ = cantidad de elementos no nulos en matriz rala $A$; $N$ = número de muestras en Ciencia de Datos; $b$ = tamaño de minilote; $\kappa = \lambda_{\max}/\lambda_{\min}$)._

|**Método**|**Hipótesis Clave**|**Costo por Iteración (Flops)**|**Memoria Requerida**|**Velocidad / Factor de Convergencia**|**Paralelismo**|
|---|---|---|---|---|---|
|**Jacobi**|$a_{ii} \neq 0$, $\rho(G_J) < 1$ (ej. diag. dom.)|$\approx 2\,\text{nnz}(A)$ (rala) o $2n^2$(densa)|Matriz $A$ + $2$ vectores ($2n$)|Lineal lenta: $\rho(G_J)^k$(Típico $1 - \mathcal{O}(h^2)$)|**Alto**(componentes independientes)|
|**Gauss-Seidel**|$a_{ii} \neq 0$, $\rho(G_{GS}) < 1$ (o $A$ SPD)|$\approx 2\,\text{nnz}(A)$|Matriz $A$ + $1$ vector ($n$, _in-place_)|Lineal: $\rho(G_{GS})^k$($\approx 2\times$ más rápido que Jacobi)|**Bajo** (sustitución secuencial)|
|**SOR ($\omega_{\text{opt}}$)**|$A$ SPD (o consist. ordenada), $0 < \omega < 2$|$\approx 2\,\text{nnz}(A) + 2n$|Matriz $A$ + $1$ vector ($n$)|Lineal acelerada: $(\omega_{\text{opt}}-1)^k = (1 - \mathcal{O}(h))^k$|**Bajo** (requiere coloreo Red-Black)|
|**SSOR + Chebyshev**|$A$ SPD, $0 < \omega < 2$, cotas $[\alpha, \beta]$|$\approx 4\,\text{nnz}(A) + \mathcal{O}(n)$ (doble barrido)|Matriz $A$ + $3$ vectores ($3n$)|Lineal rápida: $\sim \left(\frac{1-\sqrt{1-\rho}}{1+\sqrt{1-\rho}}\right)^k$|**Bajo** en barridos, **Alto** en suma vectorial|
|**Descenso por Gradiente (SD)**|$A$ SPD (o $L(x)$diferenciable convexa)|$2\,\text{nnz}(A) + 4n$($1$ prod. $Av$)|Matriz $A$ + $3$ vectores ($x, r, Ar$)|Lineal (sensible a zigzag): $\left(\frac{\kappa - 1}{\kappa + 1}\right)^k$|**Alto** (dominado por producto $Av$)|
|**Gradiente Conjugado (CG)**|$A$ SPD|$2\,\text{nnz}(A) + 10n$($1$ prod. $Ap$)|Matriz $A$ + $4$ vectores ($x, r, p, Ap$)|Finita ($\le n$); lineal en $\sqrt{\kappa}$: $2\left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k$|**Alto** (producto $Ap$ y _dot products_)|
|**PCG (Precondicionado)**|$A$ SPD, $M$ SPD fácil de invertir|$2\,\text{nnz}(A) + \text{costo}(M^{-1}r) + 12n$|Matriz $A$, $M$ + $5$vectores ($5n$)|Superlineal práctica: $2\left(\frac{\sqrt{\kappa(M^{-1}A)} - 1}{\sqrt{\kappa(M^{-1}A)} + 1}\right)^k$|Depende de $M^{-1}$ (Alto en Jacobi/Multimalla)|
|**Polyak / Nesterov (NAG)**|$L(x)$ diferenciable con gradiente Lipschitz|$1$evaluación de $\nabla L(x) + 4n$|$2$ a $3$vectores de pesos ($x^{(k)}, x^{(k-1)}, y^{(k)}$)|$\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^k$ (fuert. convexa) o $\mathcal{O}(1/k^2)$(NAG)|**Alto**|
|**Mini-batch SGD / Kaczmarz**|Suma finita $L(x) = \frac{1}{N}\sum \ell_i(x)$|$\mathcal{O}(b \cdot n)$por paso ($b \ll N$, indep. de $N$)|Pesos $x$ + activaciones del minilote $b$|Semiconvergencia rápida inicial + $\mathcal{O}(1/k)$|**Muy Alto** dentro del minilote $b$(GPU)|
|**ADAM**|Suma finita diferenciable (Redes Profundas)|$1$gradiente de minilote + $\approx 8n$|Pesos $x$ + $2$ vectores de momento ($m, v$) ($3n$)|Rápida adaptación por coordenada ($\mathcal{O}(1/\sqrt{k})$)|**Muy Alto**(operaciones elemento a elemento)|
|**Levenberg-Marquardt**|Mínimos cuadrados no lineales $\frac{1}{2}|r(x)|_2^2$|Ensamblar $J_k$($\mathcal{O}(mn)$) + resolver $\mathcal{O}(n^3)$|Matriz Jacobiana $J_k$ ($m \times n$) + $J_k^T J_k$($n \times n$)|Interpola entre Lineal (GD) y Cuadrática (Newton)|**Medio** (limitado por factorización de $J_k^T J_k + \lambda I$)|
|**ISTA / FISTA / ADMM**|Problema compuesto suave + $\ell_1$(LASSO)|ISTA: $2\,\text{nnz}(A) + \mathcal{O}(n)$; ADMM: sust. $L L^T$|Vectores primales/duales ($x, z, u$) ($3n$)|ISTA: $\mathcal{O}(1/k)$; FISTA: $\mathcal{O}(1/k^2)$; ADMM: $\mathcal{O}(1/k)$|**Alto** (umbral suave $S_\alpha$ es componente a componente)|
## 4.2. Guía Concisa de Aplicaciones Prácticas por Dominio
1. **Discretización de Ecuaciones en Derivadas Parciales (EDPs) y Sistemas Dispersos de Gran Escala (Elementos Finitos / Diferencias Finitas):**
	 - _Problema:_ Resolver sistemas $Ax = b$ con $n \sim 10^6 - 10^9$ donde $A$ es rala ($\text{nnz}(A) = \mathcal{O}(n)$) y SPD (ej. ecuación del calor, elasticidad estructural, potencial electrostático).
	- _Método de elección:_ **Gradiente Conjugado Precondicionado (PCG)** con precondicionador **Multimalla Algebraico (AMG)** o **Cholesky Incompleto ($\text{IC}(0)$)**. Los métodos de Jacobi y Gauss-Seidel/SSOR ya casi no se usan como resolvedores autónomos por su lentitud ($\rho \approx 1 - \mathcal{O}(h^2)$), pero son indispensables como **suavizadores (_smoothers_)** de errores de alta frecuencia dentro de cada nivel de los esquemas Multimalla. Si la matriz $A$ es rala pero **no simétrica** (ej. ecuaciones de convección-difusión o Navier-Stokes), CG se reemplaza por **GMRES**(basado en iteración de Arnoldi, Strang II.1).
2. **Tomografía Computarizada (CT), Resonancia Magnética (MRI) y Reconstrucción de Imágenes:**
	- _Problema:_ Reconstruir una imagen de píxeles $x$ a partir de proyecciones de rayos X o mediciones de frecuencia $b = Ax + \text{ruido}$, donde el sistema es masivo, mal condicionado o subdeterminado ($m < n$).
	- _Método de elección:_ El **Algoritmo de Kaczmarz** (conocido históricamente en medicina como _ART — Algebraic Reconstruction Technique_) proyecta iterativamente la imagen sobre los hiperplanos de cada rayo medido, operando fila por fila sin necesidad de cargar la matriz completa en RAM. Cuando se busca eliminar artefactos y reconstruir imágenes con menos dosis de radiación (_Compressed Sensing_), se añade regularización $\ell_1$ sobre las diferencias entre píxeles (_Total Variation_) o coeficientes de ondículas (_wavelets_), resolviéndose mediante **FISTA** o **ADMM** con el operador de umbral suave $S_\alpha$.
3. **Entrenamiento de Redes Neuronales Profundas (Deep Learning) y Grandes Modelos de Lenguaje:**
	- _Problema:_ Optimizar $d \sim 10^7 - 10^{12}$ pesos $x$ de una función altamente no lineal y no convexa $F(x, v)$ sobre millones de imágenes o textos.
	- _Método de elección:_ **Mini-batch SGD con Momento de Polyak/Nesterov** (muy utilizado en redes convolucionales de visión artificial como ResNet porque el ruido de los minilotes le permite escapar de mínimos locales agudos hacia valles planos que generalizan mejor) y **ADAM / AdamW**(indispensable en arquitecturas _Transformer_ y procesamiento de lenguaje natural, donde las capas de atención e incrustaciones tienen gradientes con escalas radicalmente dispares entre coordenadas). Al final del entrenamiento, el promediado de pesos **SWA** consolida la precisión del modelo.
4. **Robótica, Visión 3D (Bundle Adjustment) y Calibración de Parámetros Físicos:**
	- _Problema:_ Estimar parámetros geométricos o constantes cinéticas ajustando un modelo físico no lineal con moderada cantidad de variables ($n \sim 10^1 - 10^4$) pero alta exigencia de precisión numérica.
	- _Método de elección:_ **Algoritmo de Levenberg-Marquardt**. Su capacidad de transicionar automáticamente entre la seguridad global del Descenso por Gradiente lejos de la solución y la velocidad cuadrática de Gauss-Newton cerca del óptimo lo convierte en el estándar industrial para problemas de mínimos cuadrados no lineales.