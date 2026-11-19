# Componentes principales y la matriz más cercana de rango bajo
Los componentes principales de $A$ son sus vectores singulares, las columnas $u_j$ y $v_j$ de las matrices ortogonales $U$ y $V$.
>[!note] Análisis de los Componentes Principales (PCA)
>Usa los mayores $\sigma$ conectados a los primeros $u$ y $v$ para sintetizar la información en una matriz de datos. Dada una matriz $A$, extraemos la parte $A_k$ más importante (**los $\sigma$ más grandes**):
>$$A_k=\sigma_1u_1v_1^T+\cdots\sigma_ku_kv_k^T\qquad\text{de rango }(A_k)=k$$
>$A_k$ resuelve un problema de optimización de matrices. **La matriz de rango $k$ más cercana a $A$ es $A_k$**. En estadística identificamos las partes de $A$ con mayor varianza. Esto ubica a la SVD en el centro de la ciencia de datos.

>[!info] PCA es aprendizaje "no supervisado"
>El único instructor es el álgebra lineal: la SVD nos indica que elijamos $A_k$. Cuando el aprendizaje es supervisado, tenemos un gran conjunto de datos de entrenamiento. El aprendizaje profundo construye una función $F$ que clasifica correctamente la mayoría de los datos de entrenamiento. Luego aplicamos $F$ a nuevos datos, como veremos.

El análisis de los componentes principales está basado en la aproximación matricial por $A_k$. Que $A_k$ es la mejor opción de aproximación fue probado por Schmidt con su teorema desarrollado para operadores $A$ en espacios de funciones, y extendido directamente a matrices en espacios vectoriales.
>[!note] Propiedad clave de la matriz de rango $k$ especial $A_k$
>Si $B$ tiene rango $k$ entonces $||A-B||\geq||A-A_k||$.

Tres opciones para la norma de matrices $||A||$ tienen especial importancia y sus propios nombres:
- **Norma espectral**: $||A||_2=\max\frac{||Ax||}{||x||}=\sigma_1$ (a menudo llamada **norma $\mathcal{l}^2$**).
- **Norma de Frobenius**: $||A||_F=\sqrt{\sigma_1^2+\cdots+\sigma_r^2}$
- **Norma nuclear**: $||A||_N=\sigma_1+\sigma_2+\cdots+\sigma_r$ (la traza).
Estas normas tienen distintos valores aún para la matriz identidad de $n$ por $n$:
$$||I||_2=1\qquad||I||_F=\sqrt{n}\qquad||I||_N=n$$
Reemplazando $I$ por cualquier matriz ortogonal $Q$ las normas se mantienen iguales (porque todo $\sigma_i=1$):
$$||Q||_2=1\qquad||Q||_F=\sqrt{n}\qquad||Q||_N=n$$
Además, la norma espectral, la de Frobenius y la nuclear de cualquier matriz se mantienen iguales cuando $A$ se multiplica (sin importar el lado) por una matriz ortogonal.
Los valores singulares no cambian cuando $U$ y $V$ cambian a $Q_1U$ y $Q_2V$. Para matrices complejas, la palabra *unitaria* reemplaza a *ortogonal*. Entonces $\overline{Q}^TQ=I$. Estas tres normas son **unitariamente invariantes**: $||Q_1A\overline{Q}_2^T||=||A||$. La prueba de Mirsky del teorema de Eckart-Young se aplica a todas las normas unitariamente invariantes: $||A||$ se calcula a partir de $\Sigma$.
>[!info] Las normas espectral, de Frobenius y nuclear verifican la propiedad de Eckart-Young
>Las tres normas verifican $||Q_1AQ_2^T||=||A||$ para $Q_1$ y $Q_2$ ortogonales.
## Teorema de Eckart-Young: la mejor aproximación por $A_k$
### Ejemplo introductorio
Antes de enfocarnos en la prueba es útil ver lo que el teorema nos dice. En este ejemplo $A$ es diagonal y $k=2$:
La matriz de rango dos más cercana a $A=\begin{bmatrix}4 & 0 & 0 & 0\\ 0 & 3 & 0 & 0\\ 0 & 0 & 2 & 0 \\ 0 & 0 & 0 & 1\end{bmatrix}$ es $A_2=\begin{bmatrix}4 & 0 & 0 & 0\\ 0 & 3 & 0 & 0\\ 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0\end{bmatrix}$
Se puede decir que esta matriz diagonal es demasiado simpe y no típica. Pero la norma $L^2$ y la norma de Frobenius no cambian cuando la matriz $A$ se convierte en $Q_1AQ_2$ (para cualesquiera $Q_1$ y $Q_2$ ortogonales). Así, este ejemplo incluye cualquier matriz de 4 por 4 con valores singulares 4, 3, 2, 1. El teorema de Eckart-Young no sugiere mantener 4 y 3 porque son los más grandes. El error en $L^2$ es $||A-A_2||_2=2$. La norma de Frobenius es $||A-A_2||_F=\sqrt{5}$. La parte incómoda del problema es "las matrices de rango dos". Ese conjunto no es convexo. El promedio de $A_2$ y $B_2$ (ambas de rango 2) puede fácilmente tener rango 4. ¿Es posible que esta $B_2$ esté más cerca de $A$ que $A_2$?
Aunque los errores de $A-B_2$ son solo 0.5 en la diagonal principal donde $A-A_2$ tiene errores 2 y 1, $A_2$ es la mejor con rango $k=2$.
### Eckart-Young en $L^2$
Si rango $(B)\leq k$ entonces $||A-B||=\max\frac{||(A-B)x}{||x||}\geq\sigma_{k+1}$
Sabemos que $||A-A_k||=\sigma_{k+1}$. Toda la prueba de $||A-B||\geq\sigma_{k+1}$ depende de una buena elección del vector $x$ al calcular la norma $||A-B||$:
$$\text{Elija }x\neq0\ \text{de tal manera que }Bx=0\ \text{y }x=\sum_1^{k+1}c_iv_i$$
Primero, el espacio nulo de $B$ tiene dimensión $\geq n-k$, porque $B$ tiene rango $\leq k$. Segundo, las combinaciones de $v_1$ a $v_{k+1}$ producen un subespacio de dimensión $k+1$. Esos dos subespacios deben intersecarse. Como las dimensiones suman $(n-k)+(k+1)=n+1$, los subespacios deben compartir (al menos) una recta. Pensemos en dos planos que atraviesen $(0,0,0)$ en $\mathbb{R}^3$: comparten una recta ya que $2+2>3$. Elijamos un vector no nulo $x$ en esta recta.
Usemos ese $x$ para estimar la norma de $A-B$. Recuerde que $Bx=0$ y $Av_i=\sigma_iu_i$:
$$||(A-B)x||^2=||Ax||^2=\left|\left|\sum c_i\sigma_iu_i\right|\right|^2=\sum_1^{k+1}c_i^2\sigma_i^2$$
La suma es al menos tan grande como $\left(\sum c_i^2\right)\sigma_{k+1}^2$, que exactamente $||x||^2\sigma_{k+1}^2$. Esta última ecuación prueba que $||(A-B)x||\geq\sigma_{k+1}||x||$. Este vector $x$ nos da la cota inferior que queremos para $||A-B||$:
$$\frac{||(A-B)x||}{||x||}\geq\sigma_{k+1}\ \text{significa que }||A-B||\geq\sigma_{k+1}=||A-A_k||$$
### Norma de Frobenius
>[!info] Tres formas de la Norma de Frobenius
>La primera considera a $A$ como un vector largo, y toma la norma $l^2$ usual de ese vector. La segunda fórmula reconoce que en la diagonal principal de $A^TA$ se encuentra la norma $l^2$ de cada columna de $A$. La tercera forma se obtiene directamente de la SVD: la norma de Frobenius de $A=U\Sigma V^T$ no se ve afectada por $U$ y $V$, así $||A||_F^2=||\Sigma||_F^2$. Esto es $\sigma_1^2+\cdots+\sigma_r^2$.
>$$\begin{array}{l}||A||_F^2=|a_{11}|^2+|a_{12}|^2+\cdots+|a_{mn}|^2\quad\text{(todos los $a_{ij}^2$}\\ ||A||_F^2=\text{traza de }A^TA=(A^TA)_{11}+\cdots+(A^TA)_{nn}\\ ||A||_F^2=\sigma_1^2+\sigma_2^2+\dots+\sigma_r^2\end{array}$$

Por ejemplo, el elemento $1,1$ de $A^TA$ es la suma $|a_11|^2+\cdots+|a_{m1}|^2$ de la columna 1 de $A$. Así que la primer forma es lo mismo que la segunda, solo estamos tomando los números $|A_{ij}|^2$ de una columna por vez.
### Eckart-Young en la norma de Frobenius
Para la norma $||A-B||_F$ Pete Stewart encontró la siguiente prueba. Suponga que la matriz $B$ de rango $\leq k$ es la más cercana a $A$. Queremos demostrar que $B=A_k$. La sorpresa es que empezamos con la descomposición en valores singulares de $B$ (no de $A$):
$$B=U \begin{bmatrix}
	D & 0 \\0 & 0
\end{bmatrix}V^T\ \text{donde la matriz diagonal $D$ es de $k$ por $k$}$$
Esas matrices ortogonales $U$ y $V$ de $B$ no necesariamente diagonalizarán a $A$:
$$A=U\begin{bmatrix}L+E+R & F \\ G & H\end{bmatrix}V^T$$
Aquí $L$ es una matriz estrictamente triangular inferior en las primeras $k$ filas, $E$ es diagonal y $R$ es una matriz estrictamente triangular superior. Primero mostraremos que $L$, $R$ y $F$ son nulas al comparar $A$ y $B$ con esta matriz $C$ que claramente tiene rango $\leq k$:
$$C=U\begin{bmatrix}
	L+D+R & F\\
	0 & 0
\end{bmatrix}V^T$$
Esa es la idea clave de Stewart, construir $C$ con filas de ceros para mostrar su rango. Esas matrices ortogonales $U$ y $V^T$ no alteran la norma de Frobenius. Elevamos al cuadrado todos los elementos de la matriz y sumamos, teniendo en cuenta que $A-C$ tiene ceros donde $A-B$ tiene las matrices $L,R,F$:
$$||A-B||_F^2=||A-C||_F^2+||L||_F^2+||R||_F^2+||F||_F^2$$
Como $||A-B||_F^2||$ era la más pequeña posible concluimos que $L,R,F$ son cero. De manera similar encontramos que $G=0$. Hasta aquí sabemos que $U^TAV$ tiene dos bloques y $E$ es diagonal (como $D$):
$$U^TAV=\begin{bmatrix}E & 0\\ 0 & H\end{bmatrix}\qquad\text{y}\qquad U^TBV=\begin{bmatrix}D & 0\\ 0 & 0\end{bmatrix}$$

>[!note]
>Si $B$ es la matriz más cercana a $A$, entonces $U^TBV$ es la más cercana a $U^TAV$. Ahora vemos la luz. La matriz $D$ debe ser la misma que $E=\text{diag}(\sigma_1,...,\sigma_k)$.
>Los valores singulares de $H$ deben ser los $n-k$ valores singulares más pequeños de $A$. El menor error $||A-B||_F^2$ debe ser $||H||_F=\sqrt{\sigma_{k+1}^2+\cdots+\sigma_r^2}=\text{Eckart-Young}$
### Minimizar la distancia de Frobenius $||A-B||_F^2$
Podemos demostrar Eckart-Young al igualar las derivadas $||A-B||_F^2$ a cero. Toda la matriz de rango $k$ se factoriza en $B=CR=(m\times k)(k\times n)$. Por la SVD, podemos pedir que $C$ tenga $k$ columnas ortogonales (así $C^TC=\text{matriz diagonal }D$) y que $R$ tenga $k$ filas ortonormales (así $RR^T=I$). Buscamos obtener $C=U_k\Sigma_k$ y $R=V_k^T$. **Derivemos $E=||A-CR||_F^2$ para encontrar las matrices $C$ y $R$ que minimizan $E$**:
$$\frac{\partial E}{\partial C}=2(CR-A)R^T=0\qquad\left(\frac{\partial E}{\partial R}\right)^R=2(R^TC^T-A^T)C=0$$
La primer ecuación da $AR^T=CRR^T=C$ y la segunda da $R^TD=A^TC=A^TAR^T$. Como $D$ es diagonal, tenemos que **las columnas de $R^T$ son autovectores de $A^TA$. Son vectores singulares a derecha $v_j$ de $A$**.
Análogamente, las columnas de $C$ son autovectores de $AA^T:AA^TC=AR^TD=CD$. Entonces $C$ contiene los vectores singulares a izquierda $u_j$. ¿Qué vectores singulares minimizan el error $E$?
$E$ es la suma de todos los $\sigma^2$ que no estaban involucrados en $C$ y $R$. Para minimizar, esos deberían ser los menores valores singulares de $A$. Eso deja los mayores valores singulares para producir el mejor $B=CR=A_k$, con $||A-CR||_F^2=\sigma_{k+1}^2+\cdots+\sigma_r^2$.
## Análisis de componentes principales
Sea una matriz $A$ llena de datos, tenemos $n$ muestras. Para cada muestra medimos $m$ variables (p. ej., altura y peso). La matriz de datos $A_0$ tiene $n$ columnas y $m$ filas. En muchas aplicaciones esta es una matriz muy grande.
El primer paso es encontrar el promedio (la media muestral) de cada fila de $A_0$. Restamos la media a las $m$ entradas de esa fila. Ahora cada fila de la matriz centrada $A$ tiene media cero. Las columnas de $A$ son $n$ puntos de $\mathbb{R}^m$. Debido al centrado, la suma de los $n$ vectores columna es cero. Entonces el promedio de las columnas es el vector nulo.
A menudo esos $n$ puntos están agrupados cera de una recta o un plano u otro subespacio de baja dimensión de $\mathbb{R}^m$. Por ejemplo, en la siguiente imagen se muestra un conjunto típico de puntos muestrales agrupados alrededor de una recta en $\mathbb{R}^2$.
![[Pasted image 20240611204714.png]]
**Esta línea está alineada con el primer vector singular $u_1$ de $A$**, siendo el punto clave del PCA.
### La estadística detrás del PCA
Los números clave en probabilidad y estadística son la **media** y la **varianza**. La media es el promedio de los datos (en cada fila de $A_0$). Restar las media s a cada fila de $A_0$ produce una matriz $A$ centrada. Las cantidades cruciales son las "varianzas" y "covarianzas". Las varianzas son las sumas de los cuadrados de las distancias a la media en cada fila de $A$.
>[!note] Varianzas
>Elementos diagonales de la matriz $AA^T$.

Suponga que las columnas de $A$ corresponden a la edad de un niño en el eje $x$ y su altura en el eje $y$ (esas edades y alturas indican la diferencia con la edad y altura promedio). Buscamos la recta que pase lo más cerca posible de los puntos muestrales en la figura. Además, debemos tener en cuenta la distribución conjunta de edad-altura de los datos.
Esos son los productos punto (fila $i$ de $A$)$\cdot$(fila $j$ de $A$). Una covarianza alta significa que la altura aumenta con la edad (la covarianza negativa indica que una variable aumenta cuando la otra disminuye). Nuestro ejemplo tiene solo dos filas: una para la edad y otra para la altura: la matriz simétrica $AA^T$ es de 2 por 2. Al aumentar el número de $n$ niños en la muestra, dividimos por $n-1$ para que la escala de $AA^T$ sea estadísticamente correcta.
>[!note] Covarianzas
>Elementos de afuera de la diagonal de la matriz $AA^T$.
>La matriz de covarianza muestral se define como
>$$S=\frac{AA^T}{n-1}$$

El factor es $n-1$ porque un grado de libertad ya ha sido utilizado para media $=0$.
#### Ejemplo
He aquí un ejemplo con seis edades y alturas ya centrado para que cada fila sume cero:
$$A=\begin{bmatrix}
	3 & -4 & 7 & 1 & -4 & -3\\
	7 & -6 & 8 & -1 & -1 & -7
\end{bmatrix}$$
Para estos datos, la matriz de covarianza muestral se calcula fácilmente. Es definida positiva.
$$\textbf{Varianzas y covarianzas}\qquad S=\frac{1}{6-1}AA^T=\begin{bmatrix}
	20 & 25\\
	25 & 40
\end{bmatrix}$$
Los dos autovectores ortogonales de $S$ son $u_1$ y $u_2$. Son vectores singulares a izquierda (componentes principales) de $A$. El teorema de Eckart-Young establece que **el vector $u_1$ apunta haca la recta más cercana en la figura anterior. Los autovectores de $S$ son vectores singulares de $A$**.
El segundo vector singular $u_2$ será perpendicular a la recta más cercana.
>[!tip] PCA
>El PCA puede describirse usando la matriz simétrica $S=\frac{AA^T}{n-1}$ o la rectangular $A$. Sin duda $S$ es una buena matriz. Pero dados los datos en $A$, calcular $S$ serían un error computacional. Para matrices grandes, una SVD directa de $A$ es más rápida y más precisa.
### La geometría detrás del PCA
La mejor recta en el ejemplo dado
![[Pasted image 20240611211233.png]]
Resuelve un problema de **mínimos cuadrados perpendiculares**. Esto también se llama regresión ortogonal. Es diferente del método de cuadrados mínimos estándar para $n$ puntos de datos, o de la solución por cuadrados mínimos de un sistema lineal $Ax=b$. Ese clásico problema minimiza $||Ax-b||^2$. Mide distancias verticales (hacia arriba y abajo) hasta la mejor recta. Nuestro problema minimiza las distancias perpendiculares. El problema de cuadrados mínimos estándar nos lleva a un sistema lineal $A^TA\hat{x}=A^Tb$. Nuestro problema nos lleva a autovalores $\sigma^2$ y vectores singulares $u_1$ (autovectores de $S$).
**La suma de los cuadrados de las distancias de los datos a la recta $u_1$ es mínima**
Para demostrarlo separamos cada columna $a_j$ de $A$ en sus componentes en dirección $u_1$ y $u_2$:
$$\sum_1^n||a_j||^2=\sum_1^n|a_j^Tu_1|^2+\sum_1^n|a_j^Tu_2|^2$$
La suma de la izquierda del signo = está determinada por los datos. La primera suma de la derecha tiene términos $u_1^Ta_ja_j^Tu_1$, resultando en $u_1^T(AA^T)u_1$. Así cuando maximizamos esa suma en PCA eligiendo el primer autovector $u_1$ de $AA^T$, minimizamos la segunda suma. Esa segunda suma de las distancias al cuadrado desde los datos a la mejor recta (o el mejor subespacio) es la menor posible.
### El álgebra lineal detrás de PCA
El PCA es una manera de entender los datos: $n$ puntos muestrales $a_1,...,a_n$ es un espacio $m$-dimensional. El gráfico de datos está centrado: todas las filas de $A$ suman cero (A1=0$). La conexión crucial con el álgebra lineal se encuentra en los valores singulares de la matriz de covarianza muestral $S=\frac{AA^T}{n-1}$.
La varianza total de los datos viene de la norma de Frobenius (al cuadrado) de $A$:
$$\textbf{Varianza total}\quad T=\frac{||A||_F^2}{n-1}=\frac{||a_1||^2+\cdots+||a_n||^2}{n-1}$$
Esta es la traza de $S$: la suma de la diagonal. El álgebra lineal nos dice que la traza es igual a la **suma de los autovalores $\frac{\sigma_i^2}{n-1}$ de la matriz de covarianza muestral $S$**.
La traza de $S$ conecta la varianza total con la suma de las varianzas de las componentes principales $u_1,...,u_r$:
$$\textbf{Varianza total}\qquad T=\frac{\sigma_1^2+\cdots+\sigma_r^2}{n-1}$$
Exactamente como en la sumatoria de los cuadrados, la primer componente principal $u_1$ explica una fracción $\frac{\sigma_1^2}{T}$ de la varianza total. El siguiente vector singular $u_2$ de $A$ explica la siguiente fracción más grande $\frac{\sigma_2^2}{T}$. Cada vector singular hace lo que puede para capturar el significado de una matriz, y todos juntos lo logran.
>[!note] Finalidad del teorema de Eckart-Young
>El objetivo del teorema de Eckart-Young es que $k$ vectores singulares explican más sobre los datos que cualquier otro conjunto de $k$ vectores. Así está justificada nuestra elección de $u_1$ a $u_k$ como base para el espacio $k$-dimensional más cercano a los $n$ puntos muestrales.

El lector comprende que el grafico
![[Pasted image 20240611214726.png]]
muestra un agrupamiento de puntos de datos alrededor de una recta ($k=1$) en dimensión $m=2$. En problemas reales es usual que $k>1$ y $m>2$.
El "rango efectivo" de $A$ y $S$ es el número de valores singulares que están por encima del punto donde el ruido tapa la verdadera señal de los datos. Generalmente este punto es visible en un **gráfico de sedimentación** que muestra la inflexión en el decaimiento de los valores singulares $\sigma_i$ (y sus cuadrados). 
En el siguiente gráfico
![[Pasted image 20240611214939.png]]
Vemos el "codo" en el gráfico de sedimentación donde termina la señal y el ruido comienza a dominar.
En este ejemplo, el ruido proviene de los errores de redondeo en el cálculo de los valores singulares de la **Matriz de Hilbert** que es mal condicionada. El decaimiento en el resto de los verdaderos valores singulares es muy marcado. En la práctica, el ruido está en la misma matriz de datos: errores en la medición de $A_0$.
# Cuadrados Mínimos
Consideramos el problema $$Ax=b$$con $A\in\mathbb{R}^{m\times n}$. La existencia de soluciones de este problema equivale a que $b\in\left<a_1,a_2,...,a_n\right>$ ($a_i$ son las columnas de $A$) y en ese caso la unicidad depende de que el conjunto $\{a_1,a_2,...,a_n\}$ sea linealmente independiente.
En el caso en que $b\not\in\left<a_1,a_2,...,a_n\right>$ podemos intentar minimizar el residuo
$$r=b-Ax$$
Para eso nos conviene usar una norma: la **norma 2 es la mas adecuada** como veremos en breve.
El problema es entonces, hallar $x\in\mathbb{R}^n$ tal que
$$||Ax-b||_2=\min_{z\in\mathbb{R}^n}||Az-b||_2$$
o el equivalente a hallar $x\in\mathbb{R}^n$ tal que
$$||Ax-b||_2^2=\min_{z\in\mathbb{R}^n}||Az-b||_2^2$$
el ==enfoque tradicional (fuera del álgebra lineal numérica) es buscar una ecuación para el punto crítico $x$ como minimizante==.
Eso equivale a anular el gradiente de la función $F:\mathbb{R}^n\to\mathbb{R},F(z)=||Az-b||_2^2$, es decir
$$0=\frac{\partial F}{\partial z_k}|_x=\sum_{i=1}^m\frac{\partial([Az-b]_i)^2}{\partial z_k}|_x=$$
$$\sum_{i=1}^m\frac{\partial\left(\sum_{j=1}^na_{ij}z_j-b_i\right)^2}{\partial z_k}|_x=2\sum_{i=1}^m\left(\sum_{j=1}^na_{ij}x_j-b_i\right)a_{ik}$$
La validez de esta identidad para todo $1\leq k\leq n$ dice que
$$A^T(Ax-b)=0$$
Dicho de otro modo
$$A^TAx=A^Tb$$
estas se llaman las **ecuaciones normales**.
A pesar de su popularidad, esta no siempre es la mejor manera de calcular $x$.
Desde esas ecuaciones podemos deducir cuando habrá solución única de nuestro problema: si las columnas de $A$ son l.i. entonces $A^TAx=0$ implica
$$0=(A^TAx)^Tx=(Ax)^TAx=||Ax||_2^2$$
es decir $Ax=0$ y así $x=0$ por ser las columnas de $A$ l.i., por lo que el núcleo de $A^TA$ es trivial.
>[!info] Matriz de covarianza simétrica definida positiva
>En este caso $A^TA\in\mathbb{R}^{n\times n}$ resulta **simétrica definida positiva**. En particular podemos resolver el sistema usando Cholesky a un costo $\sim mn^2+\frac{1}{3}n^3$.
>En gran parte de las aplicaciones $m,n$ y las columnas de $A$ son l.i., lo cual muestra el éxito relativo de las ecuaciones normales.
# Interpolación polinomial
>[!info] Aplicación de cuadrados mínimos: aproximación polinomial
>Un caso tradicional de uso de cuadrados mínimos es en **aproximación polinomial**.
## Matrices de Vandermonde
### Introducción
>[!note] Matrices de Vandermonde
>Tipo específico de matriz con una estructura particular, donde cada fila es una progresión geométrica de los elementos de la primera columna. Estas matrices juegan un papel fundamental en el tema de la interpolación polinomial y tienen aplicaciones en varias áreas de las matemáticas y la ingeniería.
### Definición
Una matriz de Vandermonde $V$ de tamaño $n \times n$ para los puntos $x_1, x_2, \ldots, x_n$ se define como:
$$V = \begin{pmatrix}
1 & x_1 & x_1^2 & \cdots & x_1^{n-1} \\
1 & x_2 & x_2^2 & \cdots & x_2^{n-1} \\
1 & x_3 & x_3^2 & \cdots & x_3^{n-1} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & x_n & x_n^2 & \cdots & x_n^{n-1}
\end{pmatrix}$$
Cada elemento $V_{ij}$ de la matriz se puede expresar como $x_i^{j-1}$, donde $i$ es el índice de la fila y $j$ es el índice de la columna.
### Propiedades
Algunas propiedades importantes de las matrices de Vandermonde son:
1. **Determinante**: El determinante de una matriz de Vandermonde tiene una forma cerrada particular y se puede calcular como: $$\det(V) = \prod_{1 \leq i < j \leq n} (x_j - x_i)$$Este producto es cero si y solo si $x_i = x_j$ para algún $i \neq j$, lo que implica que las matrices de Vandermonde son invertibles si y solo si todos los $x_i$ son distintos.
2. **Invertibilidad**: Si los $x_i$ son distintos, la matriz de Vandermonde es invertible, lo cual es una propiedad clave en el contexto de la interpolación polinomial.
## Interpolación Polinomial
### Introducción a la Interpolación Polinomial

La interpolación polinomial es el proceso de encontrar un polinomio de grado \( n-1 \) que pase por un conjunto dado de \( n \) puntos. Dado un conjunto de puntos \((x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)\), buscamos un polinomio \( P(x) \) tal que:
$$P(x_i) = y_i \quad \text{para} \quad i = 1, 2, \ldots, n$$
### Sistema de Ecuaciones Lineales
Para encontrar los coeficientes del polinomio interpolador, se puede plantear un sistema de ecuaciones lineales. El polinomio \( P(x) \) se expresa como:
$$P(x) = a_0 + a_1 x + a_2 x^2 + \cdots + a_{n-1} x^{n-1}$$
Al evaluar este polinomio en los puntos \( x_i \), obtenemos el siguiente sistema de ecuaciones lineales:
$$\begin{cases}
a_0 + a_1 x_1 + a_2 x_1^2 + \cdots + a_{n-1} x_1^{n-1} = y_1 \\
a_0 + a_1 x_2 + a_2 x_2^2 + \cdots + a_{n-1} x_2^{n-1} = y_2 \\
\vdots \\
a_0 + a_1 x_n + a_2 x_n^2 + \cdots + a_{n-1} x_n^{n-1} = y_n \\
\end{cases}$$
### Uso de la Matriz de Vandermonde
El sistema de ecuaciones lineales anterior se puede escribir en forma matricial como:
$$V \begin{pmatrix}
a_0 \\
a_1 \\
a_2 \\
\vdots \\
a_{n-1}
\end{pmatrix} = \begin{pmatrix}
y_1 \\
y_2 \\
y_3 \\
\vdots \\
y_n
\end{pmatrix}
$$
Donde $V$ es la matriz de Vandermonde asociada a los puntos $x_i$. Si los puntos $x_i$ son distintos, la matriz $V$ es invertible, y los coeficientes del polinomio se pueden obtener resolviendo el sistema:
$$\begin{pmatrix}
a_0 \\
a_1 \\
a_2 \\
\vdots \\
a_{n-1}
\end{pmatrix} = V^{-1} \begin{pmatrix}
y_1 \\
y_2 \\
y_3 \\
\vdots \\
y_n
\end{pmatrix}$$
## Importancia en la Interpolación Polinomial
Las matrices de Vandermonde son cruciales para la interpolación polinomial debido a su estructura y propiedades:
1. **Eficiencia Computacional**: La estructura específica de la matriz de Vandermonde permite el uso de algoritmos eficientes para resolver el sistema de ecuaciones lineales asociado.
2. **Garantía de Solución**: La invertibilidad de la matriz de Vandermonde asegura que, para puntos \( x_i \) distintos, siempre existe un único polinomio interpolador.
3. **Análisis Teórico**: La forma del determinante de la matriz de Vandermonde proporciona una comprensión profunda de cuándo el sistema tiene una solución única, ayudando en el análisis teórico de la interpolación polinomial.
En resumen, las matrices de Vandermonde proporcionan una herramienta poderosa y eficiente para abordar problemas de interpolación polinomial, garantizando soluciones únicas y eficientes cuando los puntos de interpolación son distintos.
## Problema de interpolación
Dados $n+1$ pares de puntos $\{(x_0,y_0),(x_1,y_1)\cdots,(x_n,y_n)\}\in\mathbb{K}^n$, hallar un polinomio $p_n(x)$ de grado a lo sumo $n$ tal que $p_n(x_i)=y_i$.
Con los $x_1\neq x_j$ este problema siempre tiene solución única gracias a las siguientes observaciones.
Se propone $p_n(x)a_0+a_1x+\cdots+a_{n-1}x^{n-1}+a_nx^n$ y se hace $p_n(x_i)=y_i$ para todo $0\leq i\leq n$ que conduce a un sistema lineal con matriz de Vandermonde $V[x_0,x_1,\cdots,x_n]$
$$\begin{pmatrix}
	1 & x_0 & x_0^2 & \cdots & x_0^n\\
	1 & x_1 & x_1^2 & \cdots & x_1^n\\
	\vdots & \vdots & \vdots & \ddots & \vdots\\
	1 & x_n & x_n^2 & \cdots & x_n^n
\end{pmatrix}\begin{pmatrix}
	a_0\\a_1\\\vdots\\a_n
\end{pmatrix}=\begin{pmatrix}
	y_0\\y_1\\\vdots\\y_n
\end{pmatrix}$$
Se tiene
$$\det V[x_0,x_1,...,x_n]=\prod_{0\leq i<j\leq n}(x_j-x_i)$$
### Demostración por inducción
Como esto vale en $2\times2$, podemos demostrar por inducción suponiendo que vale hasta $n-1$.
Llamo $p_n(x)=\det V[x_0,x_1,...,x_{n-1},x]$, claramente es polinomio y tiene grado $n$, so coeficiente principal es (desarrollo por la primera fila y uso hipótesis inductiva):
$$\det V[x_0,x_1,...,x_{n-1}]=\prod_{0\leq i<j\leq n-1}(x_j-x_i)$$
Además, las raíces de $p_n(x)$ son exactamente $x_0,x_1,...,x_{n-1}$, luego se escribe como
$$p_n(x)=\prod_{0\leq i<j\leq n-1}(x_j-x_i)\prod_{k=0}^{n-1}(x-x_k)$$
evaluando en $x=x_n$ se obtiene el resultado mencionado.
### Enfoque mediante Lagrange
Podemos demostrar también la existencia de solución mediante [[Clase 16 - Extremos en regiones, multiplicadores en R2 y R3 (de Lagrange)#Multiplicadores de Lagrange|Lagrange]], donde:
$$L_i(x)=\frac{\prod_{0\leq j\leq n,j\neq i}(x-x_j)}{\prod_{0\leq j\leq n,j\neq i}(x_i-x_j)}$$
$$L_i(x_j)=\delta_i^j$$
$$p_n(x)=\sum_{i=0}^ny_iL_i(x)$$
La unicidad sale aparte.
## Interpolación en matrices no cuadradas
La teoría se mantiene en matrices de Vandermonde que son $(m>n)$.
$$\begin{pmatrix}
	1 & x_0 & x_0^2 & \cdots & x_0^n\\
	1 & x_1 & x_1^2 & \cdots & x_1^n\\
	\vdots & \vdots & \vdots & \cdots & \vdots\\
	\vdots & \vdots & \vdots & \cdots & \vdots\\
	\vdots & \vdots & \vdots & \cdots & \vdots\\
	1 & x_m & x_m^2 & \cdots & x_m^n
\end{pmatrix}$$
En estas matrices las columnas son l.i. porque las primeras $n+1$ filas forman una matriz de Vandermonde cuadrada invertible.
Si hacemos $n=1$
$$\begin{pmatrix}
	1 & x_0\\
	1 & x_1\\
	\vdots & \vdots\\
	1 & x_m
\end{pmatrix}\begin{pmatrix}
	a_0\\a_1
\end{pmatrix}=\begin{pmatrix}
	y_0\\y_1\\\vdots\\y_m
\end{pmatrix}$$
las ecuaciones normales conducen a la tradicional expresión
$$\begin{pmatrix}
	1 & \frac{1}{m+1}\sum_{i=0}^mx_i\\
	\frac{1}{m+1}\sum_{i=0}^mx_i & \frac{1}{m+1}\sum_{i=0}^mx_i^2
\end{pmatrix}\begin{pmatrix}
	a_0\\a_1
\end{pmatrix}=\begin{pmatrix}
	\frac{1}{m+1}\sum_{i=0}^my_i\\
	\frac{1}{m+1}\sum_{i=0}^mx_iy_i
\end{pmatrix}$$
![[Pasted image 20240608092534.png]]
## Interpretación geométrica
Intuitivamente si $A=U\Sigma V^*$ y $A\in\mathbb{R}^{m\times n}$, vimos que puede definirse la condición de $A$ como $\frac{\sigma_1}{\sigma_n}$. Si hacemos $A^TA=V\Sigma^2V^*$ vemos que en las ecuaciones normales la condición de la matriz final empeora $k(A^TA)=\left(\frac{\sigma_1}{\sigma_n}\right)^2$.
Al revisar el problema tenemos que
![[Pasted image 20240608092824.png]]
En este enfoque construimos el proyector $P$ sobre la imagen de $A$. Esto es, si $\{q_1,...,q_n\}$ es una base ortonormal de $\mathcal{S}=\text{Im}(A)$ tenemos
$$P_\mathcal{S}b=QQ^*b$$
donde $Q=(q_1|q_2|\cdots|q_n)$. La obtención de $Q\in\mathbb{R}^{m\times n}$ puede hacerse con la factorización $QR$, es decir
$$A=QR$$
con $R\in\mathbb{R}^{n\times n}$ triangular superior.
Escribimos
$$QRx=Ax=P_\mathcal{S}b=QQ^*b$$
es decir que basta resolver
$$Rx=Q^*b$$
>[!info] Complejidad algorítmica de interpolación polinomial
>Este algoritmo es más estable que el de las ecuaciones normales. Su costo está gobernado por $QR$ que es $\sim2mn^2-\frac{1}{3}n^3$.

Volvemos a las ecuaciones normales: supongamos que $A\in\mathbb{R}^{m\times n},m>n$, rango $n$ (columnas l.i.): debemos resolver
$$A^tAx=A^tb$$
si escribimos $A=U\Sigma V^*$, tenemos
$$V\Sigma^t\Sigma V^Tx=A^tAx=A^tb=V\Sigma^tU^tb$$
Como
$$\Sigma=\begin{pmatrix}
	\sigma_1 & 0 & 0 & \cdots & 0\\
	0 & \sigma_2 & 0 & \cdots & 0\\
	0 & 0 & 0 & \cdots & \sigma_n\\
	0 & 0 & 0 & \cdots & 0\\
	\vdots & \vdots & \vdots & \cdots & \vdots\\
	0 & 0 & 0 & \cdots & 0
\end{pmatrix}$$
y
$$\Sigma^t=\begin{pmatrix}
	\sigma_1 & 0 & 0 & \cdots & 0 & \cdots & 0\\
	0 & \sigma_2 & 0 & \cdots & 0 & \cdots & 0\\
	\vdots & \vdots & \vdots &\cdots & \vdots & \cdots & \vdots\\
	0 & 0 & 0 & \cdots & \sigma_n & \cdots & 0
\end{pmatrix}$$
tenemos
$$\Sigma^t\Sigma=\begin{pmatrix}
	\sigma_1^2 & 0 & 0 \cdots & 0\\
	0 & \sigma_2^2 & 0 & \cdots & 0\\
	0 & 0 & \ddots & 0 & 0\\
	0 & 0 & \ddots & 0 & \sigma_n^2
\end{pmatrix}$$
entonces
$$x=V\begin{pmatrix}
	\sigma_1^{-1} & 0 & 0 & \cdots & 0 & \cdots & 0\\
	0 & \sigma_2^{-1} & 0 & \cdots & 0 & \cdots & 0\\
	0 & 0 & 0 & \cdots & \sigma_n^{-1} & \cdots & 0
\end{pmatrix}U^tb$$
o sea
$$x=A^†b$$
y la pseudoinversa funciona como inversa.
En general, para $m>n$ podríamos tener rango menor a $n$ (digamos $r<n$).
Recordemos que de $A=U\Sigma V^t$, tenemos que para $1\leq i\leq r$
$$Av_i=U\Sigma V^tv_i=U\Sigma e_i=U\sigma_ie_i=\sigma_iU_i$$
y
$$Av_i=0$$
si $r<i\leq n$.
Luego $\text{ker}(A)=\left<v_{r+1},...,v_n\right>$ y $\text{Im}(A)=\left<u_1,...,u_r\right>$.
De hecho
$$A^†A=v\Sigma^†\Sigma v^t=v\begin{pmatrix}
	1 & 0 & \cdots & 0 & \cdots & 0\\
	0 & 1 & \cdots & 0 & \cdots & 0\\
	\vdots & \vdots & & \vdots & & \vdots\\
	0 & 0 & \cdots & 1 & \cdots & 0\\
	0 & 0 & \cdots & 0 & \cdots & 0\\
	\vdots & \vdots & & \vdots & & \vdots\\
	0 & 0 & \cdots & 0 & \cdots & 0
\end{pmatrix}$$
Queremos estudiar $A^†b$ en este caso.
El producto
$$U^tb=\begin{pmatrix}
	u_1^tb\\ u_2^tb \\ \vdots \\ u_1^tb
\end{pmatrix}$$
así
$$\Sigma^†U^tb=\begin{pmatrix}
	\sigma_1^{-1}(u_1^tb)\\
	\sigma_2^{-1}(u_2^tb)\\
	\vdots\\
	\sigma_r^{-1}(u_r^tb)\\
	0\\
	\vdots\\
	0
\end{pmatrix}\in\mathbb{R}^{n\times 1}$$
entonces
$$x=V\Sigma^†U^tb=\sum_{j=1}^{r}\sigma_j^{-1}(u_j^tb)v_j$$
Chequeamos
$$Ax=\sum_{j=1}^r\sigma_j^{-1}(u_j^tb)Av_j=\sum_{j=1}^r(u_j^tb)u_j=P_{\text{Im}(A)}b$$
y vemos que es solución del problema de cuadrados mínimos.
Si hubiera otro $\tilde{x}$ tal que $A\tilde{x}=P_{\text{Im}(A)}b$ entonces
$$A(\tilde{x}-x)=0$$
luego
$$\tilde{x}-x=\sum_{j=r+1}^n\alpha_jv_j$$
$$\tilde{x}=x+\sum_{j=r+1}^n\alpha_jv_j$$
como $x\perp v_j$ para todo $r+1\leq j\leq n$
$$||\tilde{x}||^2=\tilde{x}^t\tilde{x}=x^tx+\sum_{j=r+1}^n|\alpha_j|^2$$
o sea
$$||x||<||\tilde{x}||$$
la pseudoinversa da la solución de menor norma 2.
## Relación con regularización
En el caso de rango $r<n$ ( o muy mal condicionado) también se puede pensar el problema del modo siguiente: minimizar
$$F(z)=||Az-b||^2+\alpha^2||z||^2$$
y luego hacer $\alpha\to0$.
- El término se elige $\alpha^2$ para garantizar signo
- La expresión: $\alpha^2||z||^2$ se denomina penalización (se penalizan valores grandes de $||z||$)
- Para cada $\alpha$ tenemos una solución: la denotamos $x_\alpha$.
Decimos que $x_\alpha\to x$ es una solución dada por $A^†$. En efecto, para el mínimo de $F$ se tiene
$$\begin{array}{l}
	2(A^†(Ax_\alpha-b)+\alpha^2x_\alpha)=0\\
	(A^tA+\alpha^2I)x_\alpha=A^tb
\end{array}$$
en el caso más interesante en que $A$ tiene rango $r<n$ la matriz simétrica $A^tA$ no es invertible, pero para todo $\alpha\neq0,A^tA+\alpha^2I$ si lo es.
$$x_\alpha=(A^tA+\alpha^2I)^{-1}A^tb$$
usando SVD
$$A^tA+\alpha^2I=V\begin{pmatrix}
	\sigma_1^2+\alpha^2 & 0 & 0 & \cdots & 0 & \cdots & 0\\
	0 & \sigma_2^2+\alpha^2 & 0 & \cdots & 0 & \cdots & 0\\
	0 & 0 & \ddots & 0 & 0 & \cdots & 0\\
	0 & 0 & \cdots & \sigma_r^2+\alpha^2 & 0 & \cdots & 0\\
	0 & 0 & \cdots & 0 & \alpha^2 & \cdots & 0\\
	\vdots & \vdots & \vdots & \vdots & 0 & \ddots & \vdots\\
	0 & 0 & \cdots & 0 & 0 & \cdots & \alpha^2
\end{pmatrix}V^t$$
$$(A^tA+\alpha^2I)^{-1}=V\begin{pmatrix}
	(\sigma_1^2+\alpha^2)^{-1} & 0 & 0 & \cdots & 0 & \cdots & 0\\
	0 & (\sigma_2^2+\alpha^2)^{-1} & 0 & \cdots & 0 & \cdots & 0\\
	0 & 0 & \ddots & 0 & 0 & \cdots & 0\\
	0 & 0 & \cdots & (\sigma_r^2+\alpha^2)^{-1} & 0 & \cdots & 0\\
	0 & 0 & \cdots & 0 & \alpha^{-2} & \cdots & 0\\
	\vdots & \vdots & \vdots & \vdots & 0 & \ddots & \vdots\\
	0 & 0 & \cdots & 0 & 0 & \cdots & \alpha^{-2}
\end{pmatrix}V^t$$
$$(A^tA+\alpha^2I)^{-1}A^t=(A^tA+\alpha^2 I)^{-1}V\Sigma^tU^t=$$
$$V\begin{pmatrix}
	\frac{\sigma_1}{\sigma_1^2+\alpha^2} & 0 & 0 & \cdots & 0 & \cdots & 0\\
	0 & \frac{\sigma_2}{\sigma_2^2+\alpha^2} & 0 & \cdots & 0 & \cdots & 0\\
	0 & 0 & \ddots & 0 & 0 & \cdots & 0\\
	0 & 0 & \cdots & \frac{\sigma_r}{\sigma_r^2+\alpha^2} & 0 & \cdots & 0\\
	0 & 0 & \cdots & 0 & 0 & \cdots & 0\\
	\vdots & \vdots & \vdots & \vdots & 0 & \ddots & \vdots
\end{pmatrix}U^t$$
si hacemos $\alpha\to0$
$$(A^tA+\alpha^2I)^{-1}A^t\to A^†$$
# Integración numérica
>[!note] Integración numérica
>Proceso de producir un valor numérico para la integración de una función sobre un conjunto.

Aquí hay algunos ejemplos integrales que podrían ser calculadas por rutinas de computador adecuadas:
1. $\int_0^2e^{-x^2}dx$
2. $\int_0^1\int_0^1\sin{(xye^x)}\,dx\,dy$
3. $\int_0^1\int_0^1\tan{(xy^2)}\,dy\,dx$
4. $\int_0^\pi\cos{(3\cos{\theta})}\,d\theta$
Estos problemas de integración no son susceptibles a las técnicas aprendidas en el cálculo elemental. Esas técnicas dependen de la **antiderivación**. Por lo tanto, para encontrar el valor de $\int_a^bf(x)\,dx$ mediante cálculo, primero debemos producir una función $F$ con la propiedad de que $F'=f$. Entonces, tenemos que
$$\int_a^b f(x)\,dx = F(b) - F(a)$$
### Ejemplo 1
### Ejemplo 2
## Introducción a la interpolación
Una estrategia poderosa para calcular la integral $\int_a^b f(x)\,dx$ numéricamente es reemplazar $f$ por otra función $g$ que aproxime bien a $f$ y sea fácilmente integrable. Entonces, simplemente nos decimos a nosotros mismos que de  $f\approx g$ se sigue que
$$\int_a^b f(x)\,dx \approx \int_a^b g(x)\,dx$$
Ahora debería ocurrírsele al lector que los polinomios son buenos candidatos para la función $g$, y de hecho, $g$ puede ser un polinomio que interpola a $f$ en un cierto conjunto de nodos. Por supuesto, una aproximación polinómica de $f$ puede obtenerse de otras maneras, por ejemplo, truncando una serie de Taylor.
### Ejemplo
$$\int_0^1e^{x^2}\,dx\approx\int_0^1\sum_{k=0}^n\frac{x^{2k}}{k!}dx=\sum_{k=0}^n\frac{1}{(2k+1)k!}$$
Sin embargo, es deseable tener procedimientos generales que requieran solo evaluaciones del integrando. Los métodos basados en interpolación cumplen con este desiderátum. Las **funciones spline** (funciones definidas por partes, cuyos casos están representados por un polinomio) también pueden usarse para interpolar $f$, y las integrales de las funciones spline son fácilmente calculables.
## Integración mediante interpolación polinómica
Ahora abordamos la interpolación polinómica. Supongamos, entonces, que queremos evaluar la integral en $\int_a^b f(x)\,dx$. Podemos seleccionar nodos $x_0, x_1,..., x_n$ en $[a, b]$ y establecer un proceso de interpolación de [[Clase 16 - Extremos en regiones, multiplicadores en R2 y R3 (de Lagrange)#Multiplicadores de Lagrange|Lagrange]]. Definimos $$l_i(x) = \prod_{\substack{j=0 \\ j \neq i}}^n \frac{x - x_j}{x_i - x_j} \qquad (0 \leq i \leq n)$$Estos son los polinomios fundamentales para la interpolación. El polinomio de grado $\leq n$ que interpola $f$ en los nodos es $$p(x) = \sum_{i=0}^n f(x_i) l_i(x)$$Luego, como se mencionó anteriormente, simplemente escribimos $$\int_a^b f(x)\,dx \approx \int_a^b p(x)\,dx = \sum_{i=0}^n f(x_i) \int_a^b l_i(x)\,dx$$De esta manera, obtenemos una fórmula que puede usarse para cualquier $f$. Se lee de la siguiente manera: $$\int_a^b f(x)\,dx \approx \sum_{i=0}^n A_i f(x_i)$$donde $$A_i = \int_a^b l_i(x)\,dx$$Una fórmula de esta forma se llama **fórmula de Newton-Cotes**, si los nodos están igualmente espaciados.
## Regla Trapezoide
El caso más simple resulta si $n = 1$ y los nodos son $x_0 = a, x_1 = b$. Los polinomios fundamentales para la interpolación son entonces
$$l_0(x)=\frac{b-x}{b-a}\qquad l_1(x)=\frac{x-a}{b-a}$$
Consecuentemente, 
$$A_0 = \int_a^b l_0(x)\,dx = \frac{1}{2}(b - a) = \int_a^b l_1(x)\,dx = A_1$$
La fórmula de cuadratura correspondiente es 
$$\int_a^b f(x)\,dx \approx \frac{b - a}{2} [f(a) + f(b)]$$
Esto es conocido como la **regla del trapecio**. Es exacta para todo $f \in \Pi_1$ (es decir, polinomios de grado a lo sumo 1). Además, su término de error es
$$-\frac{1}{12}(b - a)^3 f''(\xi)$$
donde $\xi \in (a, b)$. Esto se puede determinar integrando el término de error en la aproximación polinómica $f(x) - p_1(x) = \frac{f''(\xi_x)(x - a)(x - b)}{2}$​, y empleando el Teorema del Valor Medio para Integrales.
Si el intervalo $[a,b]$ se particiona así
$$a = x_0 < x_1 < \cdots < x_n = b$$
entonces la regla del trapecio se puede aplicar a cada subintervalo. Aquí los nodos no están necesariamente espaciados uniformemente. Así, obtenemos la **regla del trapecio compuesta**,
$$\int_a^b f(x)\,dx = \sum_{i=1}^n \int_{x_{i-1}}^{x_i} f(x)\,dx \approx \frac{1}{2} \sum_{i=1}^n (x_i - x_{i-1}) [f(x_{i-1}) + f(x_i)]$$
(Una regla compuesta es una obtenida aplicando una fórmula de integración para un solo intervalo a cada subintervalo de un intervalo particionado). La regla del trapecio compuesta también surge si el integrando $f$ se reemplaza por una función spline interpolante de grado 1; es decir, una línea quebrada.
Con espaciamiento uniforme $h = \frac{b - a}{n}$ y $xi​=a+ih$, la regla del trapecio compuesta toma la forma $$\int_a^b f(x)\,dx \approx \frac{h}{2} \left[f(a) + 2 \sum_{i=1}^{n-1} f(a + ih) + f(b)\right]$$o
$$\int_a^b f(x)\,dx \approx h \sum_{i=0}^n{}'' f(a + ih)$$
donde la doble prima en la sumatoria indica que los primeros y últimos términos en la suma deben ser divididos por la mitad. El término de error para la regla del trapecio compuesta es
$$-\frac{1}{12}(b - a)h^2 f''(\xi)$$
donde $\xi \in (a, b)$. Se obtiene sumando el término de error en cada subintervalo y usando el hecho de que hay un punto $\xi \in [a, b]$ para el cual $f''(\xi) = \frac{1}{n}\sum_{i=1}^n f''(\xi_i)$ donde $\xi_i \in (x_{i-1}, x_i)$ y $\frac{1}{n} = \frac{b - a}{h}$​; es decir, el valor promedio.
### Ejemplo
Si tomamos $n=2$ y $[a,b]=[0,1]$ en el procedimiento de Newton-Cotes, obtenemos otra fórmula:
$$\int_0^1f(x)\,dx\approx\frac{1}{6}f(0)+\frac{2}{3}f(\frac{1}{2})+\frac{1}{6}f(1)$$
Esta fórmula puede obtenerse primero obteniendo los tres polinomios fundamentales para los nodos $0,\frac{1}{2},1$. Ellos son:
$$l_0(x)=2\left(x-\frac{1}{2}\right)(x-1)\qquad l_1(x)=-4x(x-1)\qquad l_2(x)=2x\left(x-\frac{1}{2}\right)$$
Entonces
$$A_0=\int_0^1l_0(x)\,dx=\frac{1}{6}$$
y así sucesivamente.
## Método de Coeficientes Indeterminados
En el proceso por el cual se derivó la fórmula $\int_a^b f(x)\,dx \approx \sum_{i=0}^n A_i f(x_i)$ (3), vemos de inmediato que es exacta para todos los polinomios de grado $\leq n$. Por otro lado, supongamos que se nos da una fórmula (3) y se nos dice solo que es exacta para todos los polinomios de grado $\leq n$. ¿Es cierta la siguiente ecuación?
$$A_i = \int_a^b l_i(x)\,dx$$
Sí, y vemos esto de hecho en que la fórmula (3) debería integrar correctamente cualquier $l_j$​. Por lo tanto,
$$\int_a^b l_j(x)\,dx = \sum_{i=0}^n A_i l_j(x_i) = A_j$$
Aquí, por supuesto, usamos dos propiedades esenciales de los polinomios fundamentales: $l_j$​ es un polinomio de grado $\leq n$ y $l_i(x_j) = \delta_{ij}$.
Las observaciones recién hechas nos permiten llegar a fórmulas como (3) mediante el método de los coeficientes indeterminados. Para ilustrar, vamos a derivar $\int_0^1 f(x)\,dx \approx \frac{1}{6} f(0) + \frac{2}{3} f(\frac{1}{2}) + \frac{1}{6} f(1)$ de esta manera. Buscamos una fórmula
$$\int_0^1 f(x)\,dx \approx A_0 f(0) + A_1 f(\frac{1}{2}) + A_2 f(1)$$
que sea exacta para todos los polinomios de grado $\leq 2$. Al usar como funciones de prueba los polinomios $f(x)=1$, $x$ y $x^2$, en orden, obtenemos
$$\begin{array}{l}
	1 = \int_0^1 dx = A_0 + A_1 + A_2\\
	\frac{1}{2} = \int_0^1 x\,dx = \frac{1}{2} A_1 + A_2\\
	\frac{1}{3} = \int_0^1 x^2\,dx = \frac{1}{4} A_1 + A_2
\end{array}$$
La solución de este sistema de tres ecuaciones simultáneas es $A_0 = \frac{1}{6}$​, $A_1 = \frac{2}{3}$​ y $A_2 = \frac{1}{6}$​. Dado que la fórmula es lineal, producirá valores exactos de integrales para cualquier polinomio cuadrático de la forma $f(x) = c_0 + c_1 x + c_2$.
## Regla Simpson
Cálculos similares para un intervalo arbitrario $[a,b]$ llevan a la conocida regla de Simpson:
$$\int_a^b f(x)\,dx\approx\frac{b-a}{6}\left[f(a)+4f\left(\frac{a+b}{2}\right)+f(b)\right]$$
Por la forma en que se derivó, sabemos que la regla de Simpson es exacta para todos los polinomios de grado $\leq 2$. Es, inesperadamente, exacta para todos los polinomios de grado $\leq 3$. El término de error asociado con la regla de Simpson es
$$-\frac{1}{90}(b-a)^5f^{(4)}(\xi)$$
para algún $\xi \in (a,b)$. Si $h=\frac{b-a}{2}$, esta fórmula de integración numérica toma la forma
$$\int_a^{a+2h}f(x)\,dx\approx\frac{h}{3}[f(a)+4f(a+h)+f(a+2h)]$$
Usando el [[Clase 13 - Polinomio de Taylor y aproximaciones|Teorema de Taylor]], podemos escribir el lado derecho como
$$2hf(a)+2h^2f'(a)+\frac{4}{3}h^3f''(a)+\frac{2}{3}h^4f'''(a)+\frac{100}{3 \cdot 5}h^5f^{(4)}(a)+\cdots$$
Ahora sea
$$F(x)=\int_a^x f(t)\,dt$$
Por el Teorema Fundamental del Cálculo, $F'=f$, y con el Teorema de Taylor podemos escribir el lado izquierdo en la regla de Simpson como $F(a+2h)$ o
$$2hf(a)+2h^2f'(a)+\frac{4}{3}h^3f''(a)+\frac{2}{3}h^4f'''(a)+\frac{32}{5!}h^5f^{(4)}(a)+\cdots$$
Combinando estas dos expansiones, tenemos
$$\int_a^{a+2h}f(x)\,dx=\frac{h}{3}[f(a)+4f(a+h)+f(a+2h)]-\frac{1}{90}h^5f^{(4)}(a)-\cdots$$
Y el término de error sigue.
### Composición de la regla de Simpson
Una **regla de Simpson compuesta** usando un número par de subintervalos es a menudo utilizada. Sea $n$ par, y definimos
$$x_i=a+ih \qquad h=\frac{b-a}{n} \qquad (0\leq i\leq n)$$
Entonces
$$\int_a^b f(x)\,dx=\int_{x_0}^{x_2} f(x)\,dx+\int_{x_2}^{x_4} f(x)\,dx+\cdots+\int_{x_{n-2}}^{x_n} f(x)\,dx$$
$$=\sum_{i=1}^{n/2}\int_{x_{2i-2}}^{x_{2i}} f(x)\,dx$$
Ahora aplicamos la regla de Simpson a cada subintervalo, resultando en la fórmula
$$\int_a^b f(x)\,dx \approx \frac{h}{3}\sum_{i=1}^{n/2} [f(x_{2i-1})+4f(x_{2i-1})+f(x_{2i})]$$
El lado derecho de esta fórmula debe ser calculado de la siguiente manera, para evitar la repetición de términos:
$$\frac{h}{3}\left[f(x_0)+2\sum_{i=2}^{n/2} f(x_{2i-2})+4\sum_{i=1}^{n/2} f(x_{2i-1})+f(x_n)\right]$$
El término de error para esta fórmula es
$$-\frac{1}{180}(b-a)h^4f^{(4)}(\xi)$$
para algún $\xi \in (a,b)$.
## Fórmula de integración general
El procedimiento que llevó a la fórmula de Newton-Cotes puede ser utilizado para producir una fórmula de integración más general del tipo:
$$\int_a^b f(x) w(x)\, dx \approx \sum_{i=0}^n A_i f(x_i)$$
donde $w$ puede ser cualquier función de peso fija. La única modificación necesaria es reemplazar $\int_a^b f(x) \, dx \approx \sum_{i=0}^n A_i f(x_i)$ por:
$$A_i = \int_a^b l_i(x) w(x) \, dx$$
Esto se ilustra en el siguiente ejemplo.
### Ejemplo
Encuentra una fórmula
$$\int_{-\pi}^\pi f(x)\cos{x}\,dx\approx A_0f\left(-\frac{3}{4}\pi\right)+A_1f\left(-\frac{1}{4}\pi\right)+A_2f\left(\frac{1}{4}\pi\right)+A_3f\left(\frac{3}{4}\pi\right)$$
que sea exacta para polinomios de grado tres.
#### Solución
Claramente, un polinomio de grado tres es una combinación lineal de los cuatro monomios $1$, $x$, $x^2$ y $x^3$. Así, podemos determinar los coeficientes sustituyendo $f(x)=x^j (0 \leq j \leq 3)$ y resolviendo las cuatro ecuaciones lineales resultantes. Por simetría, $A_0 = A_3$ y $A_1 = A_2$; de esta manera, el trabajo puede simplificarse:
$$0 = \int_{-\pi}^\pi 1 \cos{x} \, dx = 2A_0 + 2A_1$$
$$-4\pi = \int_{-\pi}^\pi x^2 \cos{x} \, dx = 2A_0 \left( \frac{3}{4}\pi \right)^2 + 2A_1 \left( \frac{1}{4}\pi \right)^2$$
La solución es $A_1 = A_2 = -A_0 = -A_3 = \frac{4}{\pi}$, y la fórmula es
$$\int_{-\pi}^\pi f(x)\cos{x}\,dx \approx \frac{4}{\pi} \left[ -f\left(-\frac{3}{4}\pi\right) + f\left(-\frac{1}{4}\pi\right) + f\left(\frac{1}{4}\pi\right) - f\left(\frac{3}{4}\pi\right) \right]$$
## Cambio de intervalos
De una fórmula para la integración numérica en un intervalo, podemos derivar una fórmula para cualquier otro intervalo haciendo un cambio lineal de variable. Si la primera fórmula es exacta para polinomios de cierto grado, lo mismo será cierto para la segunda. Veamos cómo se logra esto.
Supongamos que se nos da una fórmula para la integración numérica:
$$\int_c^d f(t)\,dt \approx \sum_{i=0}^n A_i f(t_i)$$
No nos importa de dónde provenga esta fórmula; sin embargo, asumamos que es exacta para todos los polinomios de grado $\leq m$. Si se necesita una fórmula para otro intervalo, digamos $[a, b]$, primero definimos una función lineal $\lambda$ de $t$ tal que si $t$ recorre $[c, d]$, entonces $\lambda(t)$ recorrerá $[a, b]$. La función $\lambda$ se da explícitamente por
$$\lambda(t) = \frac{b-a}{d-c}t + \frac{ad-bc}{d-c}$$
Ahora, en la integral
$$\int_a^b f(x)\,dx$$
hacemos el cambio de variable $x = \lambda(t)$. Entonces $dx = \lambda'(t)\,dt = (b-a)(d-c)^{-1}\,dt$, y así tenemos
$$\int_a^b f(x)\,dx = \frac{b-a}{d-c} \int_{\lambda^{-1}(a)=c}^{\lambda^{-1}(b)=d} f(\lambda(t))\,dt \approx \frac{b-a}{d-c} \sum_{i=0}^n A_i f(\lambda(t_i))$$
Por lo tanto, tenemos
$$\int_a^b f(x)\,dx \approx \frac{b-a}{d-c} \sum_{i=0}^n A_i f\left(\frac{b-a}{d-c}t_i + \frac{ad-bc}{d-c}\right)$$
Observa que, debido a que $\lambda$ es lineal, $f(\lambda(t))$ es un polinomio en $t$ si $f$ es un polinomio, y los grados son los mismos. Por lo tanto, la nueva fórmula es exacta para polinomios de grado $m$. Para la regla de Simpson, este procedimiento produce
$$\int_a^b f(x)\,dx \approx \frac{b-a}{6}\left[f(a) + 4f\left(\frac{a+b}{2}\right) + f(b)\right]$$
como consecuencia de
$$\int_0^1 f(x)\,dx \approx \frac{1}{6}f(0) + \frac{2}{3}f\left(\frac{1}{2}\right) + \frac{1}{6}f(1)$$
usando $\lambda(t) = (b-a)t + a$.
## Análisis de error
Para evaluar el error involucrado en la integración numérica $\int_a^bf(x)\,dx\approx\sum_{i=0}^nA_if(x_i)$, utilizamos el término de error en la interpolación polinómica. Recordemos, $f(x)-p(x)=\frac{1}{(n+1)!}f^{(n+1)}(\xi_x)\prod_{i=0}^n(x-x_i)$, que si $p$ es el polinomio de grado $\leq n$ que interpola $f$ en $x_0,x_1,...,x_n$, y si $f^{(n+1)}$ es continua, entonces
$$f(x)-p(x)=\frac{1}{(n+1)!}f^{(n+1)}(\xi_x)\prod_{i=0}^n(x-x_i)$$
Consecuentemente,
$$\int_a^bf(x)\,dx-\sum_{i=0}^nA_if(x_i)=\frac{1}{(n+1)!}\int_a^bf^{(n+1)}(\xi_x)\prod_{i=0}^n(x-x_i)\,dx$$
Si $|f^{(n+1)}(x)|\leq M$ en $[a,b]$, entonces
$$\left|\int_a^bf(x)\,dx-\sum_{i=0}^nA_if(x_i)\right|\leq\frac{M}{(n+1)!}\int_a^b\prod_{i=0}^n|x-x_i|\,dx$$
La elección de nodos que hace que el lado derecho de esta desigualdad sea lo más pequeño posible se conoce como
$$x_i=\frac{a+b}{2}+\frac{b-a}{2}\cos{\left[\frac{(i+1)\pi}{n+2}\right]}\qquad(0\leq i\leq n)$$
Si el intervalo es $[-1,1]$, estos nodos tienen la forma más simple
$$x_i=\cos{\left[\frac{(i+1)\pi}{n+2}\right]}\qquad(0\leq i\leq n)$$
Estos son los ceros de la función
$$U_{n+1}(x)=\frac{\sin{[(n+2)\theta]}}{\sin{\theta}}\qquad(x=\cos{\theta})$$
La función $U_{n+1}$ se conoce como un **polinomio de Chebyshev de segunda clase**. Como se define, $U_{n+1}(x)$ no es un polinomio mónico; es decir, su coeficiente principal no es 1. De hecho, el coeficiente de $x^{n+1}$ en $U_{n+1}$ es $2^{n+1}$. Así que
$$2^{-(n+1)}U_{n+1}(x)=(x-x_0)(x-x_1)\cdots(x-x_n)$$
con $x_i=\cos{[(i+1)\pi/(n+2)]}$. Una computación luego muestra que
$$\int_{-1}^1|(x-x_0)(x-x_1)\cdots(x-x_n)|\,dx=2^{-n}$$
Por lo tanto, con estos nodos y los coeficientes adecuados $A_i$ tenemos
$$\left|\int_{-1}^1f(x)\,dx-\sum_{i=0}^nA_if(x_i)\right|\leq\frac{M}{(n+1)!2^n}$$
### Teorema de Minimización de Polinomios Mónicos
>[!note] Teorema de Minimización de Polinomios Mónicos
>Entre los polinomios mónicos $p$ de grado $n$, aquel para el cual $\int_{-1}^1|p(x)|,dx$ es menor es $2^{-n}U_n$.
#### Demostración
Primero, probamos la relación de ortogonalidad
$$\int_{-1}^1U_m(x)\text{sgn}[U_n(x)]\,dx=0\qquad(0\leq m<n)$$
Sea $I$ esta integral. Hacemos el cambio de variable $\cos{\theta}=x$ para obtener
$$I=\int_0^\pi\frac{\sin(m+1)\theta}{\sin{\theta}}\text{sgn}\left[\frac{\sin{(n+1)}\theta}{\sin{\theta}}\right]\sin{\theta}\,d\theta$$
$$\sum_{k=0}^n(-1)^k\int_{k\phi}^{(k+1)\phi}\sin{(m+1)}\theta\,d\theta\qquad\phi=\frac{\pi}{n+1}$$
$$=(m+1)^{-1}\sum_{k=0}^n(-1)^{k+1}[\cos{(m+1)}(k+1)\phi-\cos{(m+1)}k\phi]$$
Sea $\alpha=(m+1)\phi+\pi$. Entonces (usando $\mathcal{Re}$ para significar parte real) tenemos
$$(m+1)I-\sum_{k=0}^n[\cos{(k+1)\alpha+\cos{k}\alpha]}=\mathcal{Re}\left\{\sum_{k=0}[e^{i\alpha(k+1)}+e^{i\alpha k}]\right\}$$
$$=\mathcal{Re}\left\{\frac{e^{i\alpha(n+2)}-e^{i\alpha}+e^{i\alpha n}-1}{e^{i\alpha}-1}\right\}$$
$$=\frac{\mathcal{Re}\{(e^{-i\alpha}-1)(e^{i\alpha(n+2)}-e^{i\alpha}+e^{i\alpha n}-1)\}}{|e^{i\alpha}-1|^2}$$
Aquí, el numerador es
$$\mathcal{Re}\{e^{i\alpha n}-e^{i\alpha(n+2)}\}=\cos{n\alpha}-\cos{(n+2)\alpha}$$
$$=\cos{[(n+1)\alpha-\alpha]}-\cos{[(n+1)\alpha+\alpha]}$$
$$=\cos{(k\pi-\alpha)}-\cos{(k\pi+\alpha)}$$
donde $k=m+n+2$. Para completar la prueba, sea $p$ un polinomio mónico de grado $n$. Puede expresarse como
$$p=2^{-n}U_n+a_{n-1}U_{n-1}+\cdots+a_0U_0$$
Por lo tanto, por la relación de ortogonalidad,
$$\int_{-1}^1|p|\,dx\geq\int_{-1}^1p\,\text{sgn}[U_n]\,dx=2^{-n}\int_{-1}^1U_n\text{sgn}[U_n]\,dx=2^{-n}\int_{-1}^1|U_n|\,dx$$