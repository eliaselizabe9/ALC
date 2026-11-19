# Descomposición de Schur (5.2 kincaid)
## Similitud de autovalores en matrices similares
>[!note] Teorema de similitud de autovalores para matrices similares
>Sean $A$ y $B$ matrices similares tales que
>$$B=PAP^{-1}$$
>Entonces, ambos poseen autovalores similares.
### Demostración
Podemos ver que $A$ y $B$ tienen el mismo polinomio característico. Entonces,
$$\begin{aligned}
	\det(B-\lambda I)&=\det(PAP^{-1}-\lambda I)\\
	&=\det[P(A-\lambda I)P^{-1}]\\
	&=\det P\cdot\det(A-\lambda I)\cdot\det P^{-1}\\
	&=\det(A-\lambda I)
\end{aligned}$$
Esto se obtiene aplicando dos propiedades básicas:
- El determinante de un producto de dos matrices es igual al producto de sus determinantes.
- El determinante de la inversa de una matriz es el inverso multiplicativo del determinante de la matriz.
## Teorema de Schur y Descomposición de Schur
A partir del teorema anterior, tenemos la sugerencia de que una posible estrategia para hallar los autovalores de una matriz $A$ es mediante su transformación a una matriz $B$. Si $B$ es más simple que $A$, el cálculo de sus autovalores puede ser más sencillo. Específicamente, suponiendo que $B$ es triangular, entonces tenemos que sus autovalores son simplemente los elementos de la diagonal (y por lo tanto, obtenemos los de $A$).
El Teorema de Schur utiliza una matriz $U$ unitaria tal que $UU^*=I$, análogo a $PP^{-1}$. Este teorema establece lo siguiente:
>[!note] Teorema de Schur
>Toda matriz cuadrada es unitariamente similar a una matriz triangular.
### Demostración
Para demostrar que toda matriz cuadrada es unitariamente similar a una matriz triangular, utilizaremos el Teorema de Schur, que establece lo siguiente:
**Teorema de Schur**: Sea $A$ una matriz cuadrada de $n \times n$ con entradas complejas. Entonces, existe una matriz unitaria $U$ y una matriz triangular superior $T$ tales que:
$$A = U T U^*$$
donde $U^*$ es la transpuesta conjugada de $U$ (también llamada adjunta de $U$).
1. **Existencia de valores propios**: Dado que $A$ es una matriz compleja cuadrada, sus valores propios existen y se pueden encontrar en el campo complejo (por el teorema fundamental del álgebra). Llamemos $\lambda_1, \lambda_2, \dots, \lambda_n$ a los valores propios de $A$.
2. **Reducción de $A$ a una forma triangular**: Sabemos que, mediante el proceso de ortogonalización de Gram-Schmidt, es posible construir una base ortonormal de vectores propios generalizados para $A$. Esta base ortonormal nos permite construir la matriz unitaria $U$, de manera que:$$U^* A U = T$$donde $T$ es una matriz triangular superior. En términos de notación, esto significa que los elementos por debajo de la diagonal principal de $T$ son cero.
3. **Propiedades de matrices unitarias**: Dado que $U$ es unitaria, se cumple que $U^* U = I$, es decir, $U$ conserva la norma de los vectores y las propiedades espectrales de $A$ (sus valores propios) en la transformación. Además, la forma triangular de $T$ garantiza que sus valores propios coinciden con los de $A$.
4. **Conclusión**: La transformación $A = UTU^*$ demuestra que $A$ es unitariamente similar a una matriz triangular superior $T$, cumpliendo el enunciado del teorema.
## Descomposición de Schur
>[!note] Descomposición de Schur
>Dada una matriz $A\in\mathbb{K}^{n\times n}$ con autovalores $\lambda_1,...,\lambda_n$ (no necesariamente distintos), existen matrices $U$ unitaria y $T$ triangular superior tales que $$A=UTU^*$$y $t_{ii}=\lambda_i$. Es decir, $A$ es unitariamente equivalente a una matriz $T$ triangular superior con los autovalores de $A$ en la diagonal en cualquier orden arbitrario. Si $A\in\mathbb{R}^{n\times n}$ y los autovalores de $A$ son reales, la matriz $U$ se puede tomar real y ortogonal.

Recordemos que para una matriz $U$ unitaria, vale $U^{-1}=U^*$, por lo tanto, las matrices $A$ y $T$ son semejantes.
### Demostración
Construimos $U$ y $T$ paso a paso.
Tomamos $\lambda_1$ autovalor de $A$ y $w_1$ autovector de $A$ con autovalor $\lambda_1$ y $||w_1||_2=1$.
Completando $\{w_1\}$ a una base de $\mathbb{K}^n$ y aplicando [[Teórica 3 - Sistemas lineales y factorización#Método de ortonormalización de Gram-Schmidt|ortonormalización de Gram-Schmidt]], podemos obtener una base ortonormal de $\mathbb{C}^n$,
$$\{w^{(1)},z^{(2)},...,z^{(n)}\}$$
Construimos la matriz $U_1$ tomando estos vectores como columnas de la matriz. Como la primera columna de $AU_1$ es $\lambda_1w^{(1)}$, obtenemos
$$U_1^*AU_1=\left[\begin{array}{l|l}\lambda_1 & * \\ \hline 0 & A_1\end{array}\right]$$
Como $U_1^*AU$ y $A$ son semejantes, tienen el mismo polinomio característico. Luego $A_1\in\mathbb{C}^{(n-1)\times(n-1)}$ tiene autovalores $\lambda_2,...,\lambda_n$. Tomamos ahora $w^{(2)}\in\mathbb{C}^{n-1}$ autovector normalizado de $A_1$ correspondiente al autovalor $\lambda_2$ y repetimos el procedimiento. Construimos $U_2\in\mathbb{C}^{(n-1)\times(n-1)}$ unitaria tal que
$$U_2^*A_1U_2=\left[\begin{array}{l|l}\lambda_2 & * \\ \hline 0 & A_2\end{array}\right]$$
y definimos
$$V_2=\left(\begin{array}{l|l}1 & * \\ \hline 0 & U_2\end{array}\right)$$
Las matrices $V_2$ y $U_1V_2$ son unitarias, y $V_2^*U_1^*AU_1V_2$ tiene la forma
$$V_2^*U_1^*AU_1V_2=\left(\begin{array}{l|l}
	\begin{matrix}
		\lambda_1 & * \\ 0 & \lambda_2
	\end{matrix} & * \\
	\hline
	0 & A_2
\end{array}\right)$$
Repitiendo este procedimiento, obtenemos matrices unitarias $U_i\in\mathbb{C}^{(n-i+1)\times(n-i+1)},i=1,...,n-1$, y matrices unitarias $V_i\in\mathbb{C}^{n\times n},i=2,...,n-1$. La matriz
$$U=U_1V_2V_3...V_{n-1}$$
es unitaria y $U^*AU$ es la factorización que buscamos.
Si todos los autovalores de $A$ son reales, los autovectores pueden elegirse también reales y todas las operaciones se pueden hacer sobre los reales, lo que prueba la última información.F
# Autovalores y autovectores de matrices especiales
## Matriz ortogonal
>[!note] Matrices ortogonales
>Si $A\in\mathbb{K}^{n\times n}$ es unitaria, todos sus autovalores cumplen $|\lambda|=1$.
## Matriz normal
>[!note] Matrices normales
>Una matriz $A$ es normal si conmuta con la matriz adjunta, es decir, si $AA^*=A^*A$.
## Matrices unitariamente equivalentes
Recordemos que dos matrices $A$ y $B$ son unitariamente equivalentes si existe $U$ unitaria tal que
$$U^*AU=B$$
## Preservación de la normalidad bajo equivalencia unitaria
>[!note] Preservación de la normalidad bajo equivalencia unitaria
>Si $A$ es normal y $B$ es unitariamente equivalente a $A$, entonces $B$ es normal.
### Demostración
Tenemos que existe $U$ unitaria tal que $U^*AU=B$. Por lo tanto,
$$BB^*=U^*AUU^*A^*U=U^*AA^*U=U^*A^*UU^*AU=B^*B$$
## Matriz unitariamente diagonalizable
Decimos que una matriz $A$ es **unitariamente diagonalizable** si existe una matriz $U$ unitaria tal que $UAU^*$ es diagonal. Esto es equivalente a decir que $A$ tiene una base ortonormal de autovectores.
>[!note] Matriz unitariamente diagonalizable
>Si $A$ es normal, entonces $A$ es unitariamente diagonalizable. En particular, autovectores de $A$ correspondientes a autovalores distintos son ortogonales.
### Demostración
Por la descomposición de Schur, sabemos que existe $T$ triangular superior unitariamente equivalente a $A$. Por lo tanto, podemos probar el resultado para $T$.
Hacemos la demostración probando mediante cálculos explícitos que una matriz $T$ triangular y normal debe ser diagonal. Tenemos
$$\begin{pmatrix}
	\overline{t}_{11} & 0 & \cdots & 0 \\
	\overline{t}_{21} & \overline{t}_{22} & \cdots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	\overline{t}_{n1} & \overline{t}_{n2} & \cdots & \overline{t}_{nn}
\end{pmatrix}\begin{pmatrix}
	t_{11} & t_{12} & \cdots & t_{1n}\\
	0 & t_{22} & \cdots & t_{2n}\\
	\vdots & \vdots & \ddots & \vdots\\
	0 & 0 & \cdots & t_{nn}
\end{pmatrix}=\begin{pmatrix}
	t_{11} & t_{12} & \cdots & t_{1n}\\
	0 & t_{22} & \cdots & t_{2n}\\
	\vdots & \vdots & \ddots & \vdots\\
	0 & 0 & \cdots & t_{nn}
\end{pmatrix}\begin{pmatrix}
	\overline{t}_{11} & 0 & \cdots & 0 \\
	\overline{t}_{21} & \overline{t}_{22} & \cdots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	\overline{t}_{n1} & \overline{t}_{n2} & \cdots & \overline{t}_{nn}
\end{pmatrix}$$
Comparando la casilla $(1,1)$ del producto a ambos lados,
$$\overline{t}_{11}t_{11}=t_{11}\overline{t}_{11}+\sum_{j=2}^nt_{1j}\overline{t}_{j1}$$
Recordando que para cualquier $a\in\mathbb{C},a\overline{a}=\overline{a}a=|a|^2$, obtenemos
$$|t_{11}|^2=|t_{11}|^2+\sum_{j=2}^n|t_{1j}|^2$$
y por lo tanto
$$0=\sum_{j=2}^n|t_{1j}|^2$$
Una suma de términos no-negativos solo puede dar 0 si todos sus términos son 0, luego
$$t_{1j}=0,j=2,...,n$$
Ahora usamos que las casillas $(2,2)$ de $T^*T$ y $TT^*$ son iguales y obtenemos
$$\overline{t}_{22}t_{22}=t_{22}\overline{t}_{22}+\sum_{j=3}^nt_{2j}\overline{t}_{j2}$$
de donde concluimos
$$t_{2j}=0,j=3,...,n$$
Repitiendo este mismo argumento, para todas las casillas de la diagonal, obtenemos
$$t_{ij}=0,\forall j>i$$
y como
$$t_{ij}=0,\forall j<i$$
porque $T$ es triangular superior, concluimos que $T$ es diagonal.
## Matrices simétricas y hermitianas
>[!note] Matrices simétricas y hermitianas
>Una matriz $A\in\mathbb{R}^{n\times n}$ es simétrica $A=A^T$. Más generalmente, $A\in\mathbb{C}^{n\times n}$ es Hermitiana si $A=A^*$.
### Propiedades de las matrices hermitianas
Si $A$ es Hermitiana, entonces
#### Realidad de la forma cuadrática
$$x^*Ax\ \text{es real para todo}\ x\in\mathbb{C}^n$$
##### Demostración
Calculamos
$$\overline{x^*Ax}=(x^*Ax)^*=x^*A^*x=x^*Ax$$
Como el conjugado de $x^*Ax$ es igual a si mismo, es un número real.
#### Realidad de los autovalores
$$\text{Todos los autovalores de $A$ son reales}$$
##### Demostración
Si $Ax=\lambda x$ y suponemos $x$ normalizado $||x||_2=x^*x=1$, entonces
$$\lambda=\lambda x^*x=x^*\lambda x=x^*Ax$$
que es número real por lo visto en a primer propiedad.
#### Hermiticidad preservada bajo transformaciones unitarias
$$S^*AS\ \text{es Hermitiana para toda}\ S\in\mathbb{C}^{n\times n}$$
##### Demostración
$(S^*AS)^*=S^*A^*S=S^*AS$, y por lo tanto $S^*AS$ es Hermitiana.
### Ejemplos
Buscamos una base ortogonal de autovectores de 
$$A=\begin{pmatrix}
	-4 & 2 & -2 \\
	2 & -7 & 4 \\
	-2 & 4 & -7
\end{pmatrix}$$
Calculamos $\text{det}(A-\lambda I_3)=-(\lambda^3+18\lambda^2+81\lambda+108)=-(\lambda+12)(\lambda+3)^2$, por lo tanto hay dos autovalores $\lambda_1=-12$ y $\lambda_2=-3$.
Calculando los núcleos de $A-\lambda_i I_3,i=1,2$, obtenemos los espacios de autovectores
$$\begin{array}{}
	E_{\lambda_1}=\left<(1,-2,2)\right>\\
	E_{\lambda_2}=\left<(2,1,0),(2,0,-1)\right>
\end{array}$$
Por la [[Teórica 5 - Matrices hermitianas. Valores singulares#Preservación de la normalidad bajo equivalencia unitaria|preservación de la normalidad bajo equivalencia unitaria]], los espacios de autovalores correspondientes a distintos autovectores son ortogonales. Por lo tanto, para obtener una base de autovectores ortogonales, solo debemos ortogonalizar cada uno de los autoespacios por separado.
El autovector $(1,-2,2)$ los normalizados a $\frac{\sqrt{5}}{5}(1,-2,2)$ y para encontrar una base ortonormal de $E_{\lambda_2}$ aplicamos la [[Teórica 3 - Sistemas lineales y factorización#Método de ortonormalización de Gram-Schmidt|ortonormalización de Gram-Schmidt]] y obtenemos $E_{\lambda_2}=\left<\frac{\sqrt{3}}{3}(2,1,0),\frac{\sqrt{5}}{15}(2,-4,-5)\right>$.
La base ortogonal de autovectores es $\mathcal{B}=\{\frac{\sqrt{5}}{5}(1,-2,2),\frac{\sqrt{3}}{3}(2,1,0),\frac{\sqrt{5}}{15}(2,-4,-5)\}$.
## Matrices definidas positivas
>[!note] Matriz definida positiva
>Cuando $A$ es simétrica, con autovalores $\lambda_1,...,\lambda_n$, definimos:
>- $A$ es *definida positiva* si $\lambda_i(A)>0$ para todo $1\leq i\leq n$.
>- $A$ es *semidefinida positiva* si $\lambda_i(A)\geq0$ para todo $1\leq i\leq n$.
### Criterio de definición de una matriz real como semidefinida positiva
>[!note] Criterio de definición de una matriz real como semidefinida positiva
>Sea $A\in\mathbb{R}^{n\times n}$ es **semidefinida positiva** si y solo si $v^TAv\geq0$ para todo $v\in\mathbb{R}^n$.
#### Demostración
Si $A\begin{pmatrix}\lambda_1 & 0 \\ 0 & \lambda_2\end{pmatrix}$, y $v=(v_1,v_2)$,
$$\begin{array}{lll}
	v^TAv & = & \begin{pmatrix}v_1 & v_2\end{pmatrix}\begin{pmatrix}\lambda_1 & 0 \\ 0 & \lambda_2\end{pmatrix}\begin{pmatrix}v_1 \\ v_2\end{pmatrix} \\
	 & = & \begin{pmatrix}v_1 & v_2\end{pmatrix}\begin{pmatrix}\lambda_1v_1 \\\lambda_2v_2\end{pmatrix} \\
	 & = & (\lambda_1v_1^2+\lambda_2v_2^2)
\end{array}$$
y $\lambda_1v_1^2+\lambda_2v_2^2\geq0$ si $\lambda_1,\lambda_2\geq0$.
### Simetría y Semidefinitud Positiva de los Productos de Gram
>[!note] Simetría y Semidefinitud Positiva de los Productos de Gram
>$A^TA$ y $AA^T$ son matrices simétricas y semidefinidas positivas. 
# Valores singulares, pseudo-inversa y aplicaciones (5.4 kincaid)
## Motivación: número de condición y autovalores
>[!tip] Referencia de mal condicionamiento de una matriz
>Una matriz está mal condicionada si está cerca de una matriz singular. O equivalentemente, está cerca de una matriz que no tiene el máximo rango posible. Podemos extender esta última idea a matrices rectangulares.
>La matriz $$A=\begin{pmatrix}1 & 1 \\ 1 & 1 \\ 1 & 1+\epsilon \end{pmatrix}$$tiene rango 2 si $\epsilon\neq0$ pero está muy cerca de una matriz de rango 1 si $\epsilon$ es chico. Podemos pensar intuitivamente que la matriz está "mal condicionada", aunque la definición que vimos no se aplica para matrices rectangulares.
### Análisis de condicionamiento de una matriz
>[!info] Análisis de condicionamiento de una matriz
>Para saber si una matriz cuadrada está mal condicionada no alcanza con ver que hay un autovalor cercano a 0 o el determinante es chico.
### Ejemplo
La matriz $A=\begin{pmatrix}\epsilon & 0 & 0 \\ 0 & \epsilon & 0 \\ 0 & 0 & \epsilon \end{pmatrix}$ está bien condicionada para cualquier $\epsilon\neq0$.
Relacionamos ahora el número de condición con los autovalores de una matriz. Para matrices simétricas,
$$\mathcal{k}_2(A)=||A||_2||A^{-1}||_2=|\lambda_\max(A)|\left|\frac{1}{\lambda_\min(A)}\right|=\frac{|\lambda_\max(A)|}{|\lambda_\min(A)|}$$
Vemos que lo que determina que una ==matriz está mal condicionada es la relación entre el autovalor de mayor módulo y el autovalor de menor módulo==.
## Introducción
### Ejemplo introductorio
Generemos matrices simétricas al azar y grafiquemos la imagen del círculo unitario por la transformación lineal resultante. Para conseguir figuras uniformes, reescalamos los gráficos de modo que el círculo unitario (en azul) contenga a su imagen (en rojo). El resultado obtenido es esperable. Las imágenes son elipses debido a que cada una de estas matrices posee una base $\{v_1,v_2\}$ de autovectores ortonormales. Es decir, $Av_i=\lambda_iv_i$, por lo que cualquier punto $v$ del circuito unitario puede escribirse como $v=x_1v_1+x_2v_2$ con $1=x_1^2+x_2^2$, por lo que $Av=x_1\lambda_1v_1+x_2\lambda_2v_2$ es un punto de una elipse de semiejes a lo largo de las rectas generadas por $v_1$ y $v_2$ de longitudes $|\lambda_1|$ y $|\lambda_2|$.
![[Pasted image 20240527162727.png]]
Vemos que algo similar ocurre con matrices aleatorias, como en el siguiente gráfico.
![[Pasted image 20240527162750.png]]
En este caso no podemos garantizar una base de autovectores ortonormales, sin embargo podemos estudiar las pre-imágenes de los semiejes de las elipses (en verde y rosa respectivamente) y notamos que provienen de vectores unitarios ortogonales (con el mismo color en las guras). Es decir que nuestro experimento sugiere la existencia de dos pares de vectores ortonormales $\{v_1,v_2\}$ y $\{u_1,u_2\}$ y dos escalares $\sigma_1,\sigma_2$, que podemos suponer ordenados (i.e. $\sigma_1\geq\sigma_2\geq0$) y no negativos, cambiando si hiciera falta $u_i$ por $-u_i$, tales que $Av_i=\sigma_iw_i$.
Esta observación nos permite escribir
$$AV=U\Sigma$$
donde
$$V=\begin{pmatrix}| & |\\v_1 & v_2\\ | & |\end{pmatrix}\qquad\Sigma=\begin{pmatrix}\sigma_1 & 0 \\ 0 & \sigma_2\end{pmatrix}\qquad U=\begin{pmatrix}| & |\\u_1 & u_2\\ | & |\end{pmatrix}$$
con matrices $V,U$ ortogonales. Dicho otro modo puede factorizarse $A$ como
$$A=U\Sigma V^*$$
>[!tip] La factorización en valores es una generalización de la diagonalización de matrices simétricas
>Esta factorización generaliza en algún sentido a la obtenida al diagonalizar matrices simétricas. Si bien, no toda matriz es diagonalizable en una base ortonormal $\{v_1,v_2\}$, nuestro experimento dice que algo similar puede obtenerse si nos permiten trabajar con *bases diferentes*: una $\{v_1,v_2\}$ para el espacio de partida y otra $\{u_1,u_2\}$ para el de llegada.
>Esta factorización **siempre existe y no solo para matrices cuadradas**.
### Descomposición en valores singulares (SVD)
>[!note] Descomposición en valores singulares
>Dada $A\in\mathbb{K}^{m\times n}$, la descomposición en valores singulares es una factorización de $A$ $$A=U\Sigma V^*$$con $U\in\mathbb{K}^{m\times m}$ y $V\in\mathbb{K}^{n\times n}$ matrices unitarias y $\Sigma\in\mathbb{R}^{m\times n}$. La matriz $\Sigma$ tiene la forma $$\Sigma=\begin{pmatrix} D \\ 0_{(m-n)\times n}\end{pmatrix}$$si $m\geq n$, o la forma $$\Sigma=\begin{pmatrix}D & 0_{m\times(n-m)}\end{pmatrix}$$si $n>m$. En ambos casos $D\in\mathbb{R}^{p\times p}$, donde $p=\min\{m,n\}$, es una matriz diagonal con elementos en la diagonal $d_{11}\geq d_{22}\geq\cdots\geq d_{pp}\geq0$. Los elementos no nulos de la diagonal de $D$ se denominan valores singulares de $A$, y se notan $\sigma_i(A)=d_{ii}$.
#### Ejemplo
$$\begin{pmatrix}
	3 & 2 & 2 \\
	2 & 3 & -2
\end{pmatrix}=\begin{pmatrix}
	\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2} \\
	\frac{\sqrt{2}}{2} & -\frac{\sqrt{2}}{2}
\end{pmatrix}\begin{pmatrix}
	5 & 0 & 0 \\
	0 & 3 & 0
\end{pmatrix}\begin{pmatrix}
	\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2} & 0 \\
	\frac{\sqrt{2}}{6} & -\frac{\sqrt{2}}{6} & \frac{2\sqrt{2}}{3} \\
	\frac{2}{3} & -\frac{2}{3} & -\frac{1}{3}
\end{pmatrix}$$
## Construcción
Para construir la descomposición SVD observamos las siguientes relaciones.
Dada una matriz $A\in\mathbb{K}^{m\times n}$, si $A=U\Sigma V^*$:
- $AA^*=(U\Sigma V^*)(V\Sigma^*U^*)=U(\Sigma\Sigma^*)U^*$, con $\Sigma\Sigma^*\in\mathbb{R}^{m\times m}$ diagonal.
- $A^*A=(V\Sigma^*U^*)(U\Sigma V^*)=V(\Sigma^*\Sigma)V^*$, con $\Sigma^*\Sigma\in\mathbb{R}^{n\times n}$ diagonal.
Recordar que $A^*A$ y $AA^*$ son matrices semidefinidas positivas y por lo tanto diagonalizables con autovalores no-negativos y base ortonormal de autovectores.
En base a estas relaciones, concluimos:
- Las columnas de $U$ son una base de autovectores de $AA^*\in\mathbb{K}^{m\times m}$
- Las columnas de $V$ son una base de autovectores de $A^*A\in\mathbb{K}^{n\times n}$
- Los valores singulares en $\Sigma$ son la raíz cuadrada de los autovalores no nulos de $A^*A$ o $AA^*$ ordenados de mayor a menor.
Esto nos da un **método para calcular la descomposición**.
Primero construimos la matriz $V$ tomando una base de autovectores de $A^*A$.
Luego obtenemos las filas de $U$a partir de la relación
$$A=U\Sigma V^*$$
multiplicando por $V$:
$$AV=U\Sigma$$
En el término de izquierda tenemos
$$AV=\begin{pmatrix}
	| & & |\\
	Av_1 & \cdots & Av_n\\
	| & & |
\end{pmatrix}$$
1. Si $m\leq n$,$$AV=U\Sigma=\begin{pmatrix}| & & | \\ u_1 & \cdots & u_m \\ | & & |\end{pmatrix}\left(\begin{array}{lll|l}\sigma_1\\ & \ddots & & 0\\ & & \sigma_m\end{array}\right)=\left(\begin{array}{lll|l} | & & | \\ \sigma_1u_1 & \cdots & \sigma_mu_m & 0\\ | & & |\end{array}\right)$$Por lo tanto, para $1\leq i\leq n$, las columnas de $U$ cumplen $\sigma_iu_i=Av_i$, donde $v_i$ son las columnas de $V$. Obtenemos $$u_i=\frac{1}{\sigma_i}Av_i$$para $1\leq i\leq m$ y $\sigma_i\neq0$. Si existen valores singulares nulos, las columnas correspondientes de $U$ las elegimos completando una base ortonormal $\{u_1,...,u_m\}$ de $\mathbb{K}^{m\times m}$.
2. Si $m\geq n$, $$AV=U\Sigma=\begin{pmatrix}| & & | \\ u_1 & \cdots & u_m \\ | & & |\end{pmatrix}\begin{pmatrix}\sigma_1\\ & \ddots \\ & & \sigma_n \\ \hline & 0\end{pmatrix}=\begin{pmatrix}| & & | \\ \sigma_1u_1 & \cdots & \sigma_n u_n\\ | & & |\end{pmatrix}$$Por lo tanto, para $1\leq i\leq n$, las columnas $U$ cumplen $\sigma_iu_i=Av_i$, donde $v_i$ son las columnas de $V$. Obtenemos $$u_i=\frac{1}{\sigma_i}Av_i$$para $1\leq i\leq n$ y $\sigma_i\neq0$. El resto de las columnas de $U$ las elegimos completando una base ortonormal $\{u_1,...,u_m\}$ de $\mathbb{K}^{m\times m}$.
### Teorema de existencia de una descomposición en valores singulares
>[!note] Existencia de descomposición en valores singulares
>Dada una matriz rectangular $A\in\mathbb{K}^{m\times n}$, existen matrices
>- $U\in\mathbb{K}^{m\times m}$ unitaria,
>- $\Sigma\in\mathbb{R}^{m\times n}$, con $\Sigma=(D|0)$ o $\Sigma=\left(\frac{D}{0}\right)$, $D$ diagonal y $d_{11}\geq\cdots\geq d_{ss}\geq0$ ($s=\min\{m,n\}$),
>- $V\in\mathbb{K}^{n\times n}$ unitaria,
>tales que $A=U\Sigma V^*$. Los valores no nulos $d_{11}\geq\cdots\geq d_{ss}>0$ son únicos y se denominan *valores singulares* de $A,\sigma_i=d_{ii}$.
#### Demostración
En la construcción que hicimos, falta ver que los vectores $u_i=\frac{1}{\sigma_i}Av_i$ construidos son ortogonales. Calculamos
$$\begin{array}{lll}
	u_j^*u_i & = & \left(\frac{1}{\sigma_j}v_j^*A^*\right)\left(Av_i\frac{1}{\sigma_i}\right)\\
	 & = & \frac{1}{\sigma_j}v_j^*(A^*A)v_i\frac{1}{\sigma_i}=\frac{1}{\sigma_j}(v_j^*)(V\Sigma^*\Sigma V^*)v_i\frac{1}{\sigma_i}
\end{array}$$
y usando que $V^*v_i=e_i$, obtenemos
$$u_j^*u_i=\frac{1}{\sigma_j}(e_j^*\Sigma^*)(\Sigma e_i)\frac{1}{\sigma_i}=\frac{1}{\sigma_j}(\sigma_je_j^*)(\sigma_ie_i)\frac{1}{\sigma_i}=e_j^*e_i=\begin{cases}1 & \text{si }i=j\\0 & \text{si }i\neq j\end{cases}$$
## Valores singulares y norma-2
### Teorema de la Norma Espectral y el Valor Singular Máximo
>[!note] Teorema de la Norma Espectral y el Valor Singular Máximo
>Dada una matriz rectangular $A\in\mathbb{K}^{m\times n}$,
>$$||A||_2=\sigma_1(A)$$
#### Demostración
Recordando que multiplicar por una matriz unitaria es una transformación inversible que no cambia la norma de un vector,
$$\begin{array}{lll}
	||A||_2 & = & \sup_{||x||_2=1}||Ax||_2\\
	 & = & \sup_{||x||_2=1}||U\Sigma V^*x||_2=\sup_{||x||_2=1}||\Sigma V^*x||_2\\
	 & = & \sup_{||y||_2=1}||\Sigma y||_2=\sup_{||y||_2=1}\left(\sum_{i=1}^r\sigma_i^2|y_i|^2\right)^{\frac{1}{2}}\leq\sigma_1
\end{array}$$
### Valor Singular Mínimo y Norma 2 Mínima
>[!note] Valor Singular Mínimo y Norma 2 Mínima
>Si $A\in\mathbb{K}^{m\times n}$ y $\text{rg}(A)=n$, entonces
>$$\min_{||x||_2=1}||A||_2=\sigma_n(A)$$
### Corolario: Condición de número espectral
>[!note] Condición de número espectral
>Si $A\in\mathbb{K}^{n\times n}$ inversible,
>$$\text{cond}_2(A)=\frac{\sigma_1(A)}{\sigma_n(A)}$$
>con $r=\min\{m,n\}$.
## Distancia a matrices de menor rango
Dada una matriz $A\in\mathbb{K}^{n\times n}$, buscamos una matriz singular lo más cercana posible a $A$. En otras palabras, buscamos $B\in\mathbb{K}^{n\times n}$ de norma lo más chica posible tal que $A+B$ sea singular.
1. $A+B=A(I+A^{-1}B)$
2. Si $||A^{-1}B||<1$, entonces $I+A^{-1}B$ (y por lo tanto $A+B$) es inversible.
3. Por lo tanto, si $A+B$ es singular y $A$ es inversible, debe ser $||B||\geq \frac{1}{||A^{-1}||}$.
4. Tomando norma-2, y $A=U\Sigma V^*$ la SVD, $||A^{-1}||_2=\frac{1}{\sigma_n}$.
5. Por lo tanto, si $A+B$ es singular, debe ser $||B||\geq\sigma_n(A)$.
6. Podemos tomar $B=UEV^*$, con $E=\text{diag}(0,0,...,0,-\sigma_n)$, que cumple $||B||_2=||E||_2=\sigma_n$.
### Teorema de la Matriz Singular más cercana
>[!note] Teorema de la Matriz Singular más cercana
>Dada $A\in\mathbb{K}^{n\times n}$ inversible, con $A=U\Sigma V^*$ su descomposición SVD. Definiendo $$\Sigma_{n-1}=\begin{pmatrix}\sigma_1 \\ & \ddots \\ & & \sigma_{n-1}\\ & & & 0\end{pmatrix}$$la matriz $\tilde{A}=U\Sigma_{n-1}V^*$ es la matriz singular más cercana a $A$ en norma-2.
### Teorema de la Matriz de Rango Reducido más cercana
>[!note] Teorema de la Matriz de Rango Reducido más cercana
>Sea $A\in\mathbb{K}^{m\times n}$ de rango $k>0$, con $A=U\Sigma V^*$ su descomposición SVD. La matriz $A_1\in\mathbb{K}^{m\times n}$ de rango $k_1<k$ más cercana a $A$ en la norma 2 (de Frobenius) es $A_1=U\Sigma_1 V^*$, con $$\Sigma_1=\begin{pmatrix}\sigma_1\\ & \ddots \\ & & \sigma_{k_1}\\ & & & 0\end{pmatrix}$$la matriz igual a $\sigma$ excepto que solo se usan los valores singulares $\sigma_1,...,\sigma_{k_1}$ y los demás se reemplazan por 0.
## Aplicaciones
### Reducción de dimensionalidad
Para una población de $m$ individuos analizamos $n$ variables (features), con $m>n$, y queremos ver si hay variables redundantes o si podemos reducir la cantidad de variables sin perder mucha información.
Dada una matriz $A\in\mathbb{R}^{m\times n},m>n$, el $k$-ésimo valor singular nos dice la distancia de $A$ a una matriz de rango $k-1$. Es decir, ==mide "cuanta información" vamos a perder si utilizamos $k-1$ variables en lugar de las $n$ variables originales==.
### Pseudo-inversa
A partir de la descomposición en valores singulares $A=U\Sigma V^*$, construimos la pseudo-inversa de $A$,
$$A^†=V\Sigma^†U^*$$
donde $\Sigma^†$ se obtiene a partir de $\Sigma$ transponiéndola y reemplazando los elementos no nulos $d_{ii}$ en la diagonal por $\frac{1}{d_{ii}}$.
#### Propiedades de la pseudo-inversa
Se cumplen las siguientes propiedades:
- $AA^†A=A$ y $A^†AA^†=A^†$ (pero en general no se cumple $AA^†=I_m$ o $A^†A=I_n$).
- Si $A\in\mathbb{K}^{n\times n}$ es inversible, $A^{-1}=A^†$.
### Resolución de sistemas de ecuaciones
#### Teorema de la Solución de Mínima Norma-2 para Sistemas Lineales
>[!note] Teorema de la Solución de Mínima Norma-2 para Sistemas Lineales
>Dada la ecuación $Ax=b$, con $A\in\mathbb{K}^{m\times n}$, si el sistema tiene soluciones entonces $$z=A^†b$$es una solución y si hay infinitas soluciones, $z$ es la solución de mínima norma-2.
#### Criterio de existencia de solución mediante pseudo-inversa
Utilizando la pseudo-inversa podemos obtener también un criterio para determinar si un sistema tiene solución.
>[!note] Criterio de existencia de solución mediante pseudo-inversa
>El sistema $Ax=b$ tiene solución si y solo si $AA^†b=b$.
##### Demostración
- Si $AA^†b=b$, entonces $v=A^†b$ es solución del sistema.
- Recíprocamente, si existe $v$ tal que $Av=b$, entonces:
$$b=Av=AA^†Av=AA^†b$$
#### Teorema de Solución Generalizada para Sistemas Lineales
>[!note] Teorema de Solución Generalizada para Sistemas Lineales
>Dadas $A\in\mathbb{K}^{m\times n}$ y $b\in\mathbb{K}^n$. Si $AA^†b=b$, todos los vectores de la forma $$x=A^†b+(I_n-A^†A)y$$con $y\in\mathbb{K}^n$ son soluciones del sistema $Ax=b$. Más aún, cualquier solución puede ser escrita de esa forma.
##### Demostración
Si $x$ es de esa forma, entonces
$$\begin{array}{lll}
	Ax & = & A(A^†b+(I_n-A^†A)y)\\
	 & = & (AA^†b)+(A-AA^†A)y\\
	 & = & b+(A-A)y=b
\end{array}$$
Ahora, dado un vector $z\in\mathbb{K}^n$ tal que $Az=b$, escribimos $$z=(A^†A)z+(I_n-A^†A)z=A^†b+I_n-A^†A)z$$que es de la forma buscada.
##### Observación
- Puede haber soluciones repetidas entre las soluciones obtenidas de esta forma.
- En particular, el sistema $Ax=b$ tiene solución única si y solo si $A^†A=I_n$.
### Compresión de imágenes
Veamos como comprimir una imagen aplicando reducción de dimensionalidad por SVD.
Tomamos una matriz $A\in\mathbb{R}^{m\times n}$ , con $a_{ij}\in[0, 1]$, que representa una tonalidad de gris, siendo $0 = \text{negro}$ y $1 = \text{blanco}$.
#### Ejemplo 1
Para representar la imagen de la figura, tomamos $A\in\mathbb{R}^{25\times15}$
![[Pasted image 20240527184704.png]]
¿Cuál es el rango de $A$? ¿Cuántos valores singulares no nulos tiene $A$?
##### Eliminamos columnas de $U$
La descomposición en SVD de $A\in\mathbb{R}^{25\times15}$,
$$A=U\Sigma V^*$$
tiene 3 valores singulares no nulos, $U\in\mathbb{R}^{25\times25},\Sigma\in\mathbb{R}^{25\times15},V\in\mathbb{R}^{15\times15}$. ¿Podemos descartar filas de $U$ y $V$ sin perder información?
Recordando que
$$U\Sigma=\begin{pmatrix}
	| & & |\\
	u_1 & \cdots & u_{25}\\
	| & & |
\end{pmatrix}\begin{pmatrix}
	\sigma_1\\
	 & \ddots\\
	  & & \sigma_{15}\\
	\hline & 0
\end{pmatrix}=\begin{pmatrix}
	| & & |\\
	\sigma_1u_1 & \cdots & \sigma_{15}u_{15}\\
	| & & |
\end{pmatrix}$$
vemos que podemos tomar siempre $u_i=0$ para $i>16$ sin perder información, y en nuestro caso también $u_i=0$ para $i>3$.
##### Eliminamos columnas de $V$
Análogamente, a partir de $(\Sigma V^*)^*$ tenemos que
$$V\Sigma^*=\begin{pmatrix}
	| & & |\\
	v_1 & \cdots & v_{15}\\
	| & & |
\end{pmatrix}\left(\begin{array}{lll|l}
	\sigma_1 \\
	 & \ddots & & 0\\
	 & & \sigma_{15}
\end{array}\right)=\left(\begin{array}{lll|l}
	| & & |\\
	\sigma_1v_1 & \cdots & \sigma_{15}v_{15} & 0\\
	| & & |
\end{array}\right)$$
vemos que podemos tomar $v_j=0$ para $j>3$ sin perder información.
En forma más compacta, para no guardar tantos $0's$, podemos guardar la descomposición
$${A \atop 25\times15}={\begin{pmatrix}
	| & & |\\
	u_1 & \cdots & u_3\\
	| & & |
\end{pmatrix} \atop25\times3}{\begin{pmatrix}
	\sigma_1 & 0 & 0 \\
	0 & \sigma_2 & 0 \\
	0 & 0 & \sigma_3
\end{pmatrix} \atop 3\times3}{\begin{pmatrix}
	-v_1-\\
	-v_2-\\
	-v_3-
\end{pmatrix}\atop 3\times15}$$
##### Cantidad de datos
En total, necesitamos guardar:
- $3$ columnas de $U\to75$ valores
- $3$ columnas de $V\to45$ valores
- $3$ valores singulares
Total: $75+45+3=123$ valores.
Para guardar la matriz original necesitamos $25\times15=375$ valores.
#### Ejemplo 2
Tomamos ahora la imagen de un árbol $A\in\mathbb{R}^{1082\times2000\times3}$
![[Pasted image 20240527190403.png]]
Para imágenes en color se utilizan tres matrices con tonalidades RGB (red, green, blue), pero vamos a utilizar solo tonalidades de gris, tomamos $A\in\mathbb{R}^{1080\times2000}$.
![[Pasted image 20240527190434.png]]
Calculamos la descomposición $A=U\Sigma V^*$ y graficamos los logaritmos de los valores singulares. Podemos tomar los primeros 987 valores singulares sin casi perder información, pero no es suficiente compresión.
![[Pasted image 20240527190528.png]]
Probamos tomar las primeros 200 valores singulares. Es decir, tomamos
$$\tilde{A}={\tilde{U}\atop1080\times200}\quad{\tilde{\Sigma}\atop200\times200}\quad{\tilde{V}^*\atop200\times2000}$$
y obtenemos la imagen
![[Pasted image 20240527190720.png]]