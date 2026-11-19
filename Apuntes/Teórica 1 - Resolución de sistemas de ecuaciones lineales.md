# Repaso: Vectores y matrices
## Vectores
>[!note] Vector
>Sucesión de $n$ ($n\in\mathbb{N}$) números reales o complejos $v=(v_1,...,v_n)$.

Los vectores pueden pertenecer a $\mathbb{R}^n$ o a $\mathbb{C}^n$. Pero cuando los resultados valen para ambos conjuntos, nos referimos a estos como $\mathbb{K}$.
Algunos ejemplos son:
- $v=(1-2i,2+3i)\in\mathbb{C}^2$
- $u=(2,-1,\pi,3/2)\in\mathbb{R}^4$
### Definición de vectores en Python
En **Python** definimos vectores con el comando `array` del paquete `numpy`. Por ejemplo:
``` Python
import numpy as np

v1 = np.array([10,5,-7,1])
print(v1)

v2 = np.array([5,0,3/2,2])
print(v2)
```
## Matrices
>[!note] Matrices
>Denotamos $\mathbb{K}^{m\times n}$ al **espacio de matrices** de números reales o complejos de $m$ filas y $n$ columnas, que podemos considerar como un vector de $m\times n$ coordenadas agrupadas en $m$ filas de $n$ coordenadas.

Por ejemplo:
- $\begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6\end{pmatrix}\in\mathbb{R}^{2\times3}$
- $\begin{pmatrix}1 \\ 7 \\ 0.6 \\ -1 \end{pmatrix}\in\mathbb{R}^{4\times1}$
- $\begin{pmatrix} 1 & 5+i \\ 7-2i & 0 \\ i & 3 \end{pmatrix}\in\mathbb{C}^{3\times2}$

>[!tip] Vectores como matrices columna
>Es común considerar a los vectores en $\mathbb{K}^n$ como matrices columna. Es decir, al vector $v=(1,2,3)$ lo pensamos como matriz $\begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}$.
#### Definición de matrices en Python
Al igual que con los vectores, en **Python** definimos matrices con `array`, como un arreglo de $m$ vectores, donde cada vector es una fila. Por ejemplo:
``` Python
import numpy as np
A = np.array([[1,2],[3,4]]) # Matriz de 2 x 2
print("A = \n", A)

B = np.array ([[1,2,3,4],[7,1,2,-1]]) # Matriz de 2 x 4
print("B = \n", B)

C = np.array([[1],[7],[1/3]]) # Matriz columna de 1 x 3
print("C = \n", C)
```
### Matrices canónicas
>[!note] Matrices canónicas
>Sean $n,m\in\mathbb{N}$, y sean $1\leq k\leq n,1\leq l\leq m$. Se define la matriz $E^{kl}\in K^{n\times m}$ como
>$$(E^{kl})_{ij}=\begin{cases}1 & \text{si $i=k,j=l$} \\ 0 & \text{si no}\end{cases}$$
### Matrices como espacios vectoriales
Sean $A,B\in K^{n\times m}$. Entonces $A=B$ si y sólo si $A_{ij}=B_{ij}$ para cada $1\leq i\leq n,1\leq j\leq m$. De este modo, ==podemos definir una operación (suma) en $K^{n\times m}$ y una acción de $K$ en $K^{n\times m}$ que transforman a este conjunto en un $K$-espacio vectorial.==

>[!note] Suma de matrices y producto por escalares
>$+:K^{n\times m}\times K^{n\times m}\to K^{n\times m}, (A+B)_{ij}=A_{ij}+B_{ij}\quad(1\leq i\leq n,1\leq j\leq m)$
>$\cdot:K\times K^{n\times m}\to K^{n\times m},(\lambda\cdot A)_{ij}\quad(1\leq i\leq n,1\leq j\leq m)$

Ya definidas estas operaciones, es fácil verificar que $(K^{n\times m},+,\cdot)$ es un $K$-espacio vectorial.
### Producto de matrices
>[!note] Producto de matrices
>Sean $A\in K^{n\times m}$ y $B\in K^{m\times r}$. Se define el producto de $A$ por $B$ como la matriz $C\in K^{n\times r}$ tal que
>$$C_{ij}=\sum_{k=1}^m A_{ik}B_{kj}\qquad 1\leq i\leq n,1\leq j\leq r$$
#### Propiedades
1. **Propiedad asociativa**: dadas $A\in K^{n\times m}$, $B\in K^{m\times r}$ y $C\in K^{r\times s}$, se tiene que $(A\cdot B)\cdot C=A\cdot(B\cdot C)$.
2. **Elemento neutro**: Para cada $n\in\mathbb{N}$, sea $I_n\in K^{n\times n}$ definida por $(I_n)_{ij}=\begin{cases} 1 & \text{si } i=j \\ 0 & \text{si } i\neq j \end{cases}$ Entonces, si $A\in K^{n\times m}$, se verifica que $I_n\cdot A=A\cdot I_m=A$
3. **Propiedades distributivas**:
	1. Si $A\in K^{n\times m}$ y $B,C\in K^{m\times r}$, entonces $A\cdot(B+C)=A\cdot B+A\cdot C$
	2. Si $A,B\in K^{n\times m}$ y $C\in K^{m\times r}$, entonces $(A+B)\cdot C=A\cdot C + B\cdot C$
>[!warning] Propiedades que no se cumplen para el producto de matrices
>- **El producto de matrices no es conmutativo**. Aún en el caso de matrices cuadradas, en el que siempre se pueden calcular $A\cdot B$ y $B\cdot A$, en general se tiene que $A\cdot B\neq B\cdot A$.
>- Si consideramos que $A\cdot B=0$, esto no implica que $A=0$ o $B=0$.
### Producto de matrices $K^{n\times n}$
#### Define un anillo
El producto de matrices está definido para cualquier par de matrices en $K^{n\times n}$ y, por lo tanto, se tiene una operación "producto" en $K^{n\times n}$ para cada $n\in\mathbb{N}$. Por lo tanto, $(K^{n\times n},+,\cdot)$ es un anillo.
#### Definición Fundamental de un $K$-álgebra
Sea $K$ un cuerpo y sea $A$ un conjunto con dos operaciones, $+$ y $\cdot$, y una acción $\cdot_K$ de $K$ en $A$ tales que:
1. $(A,+,\cdot)$ es un anillo
2. $(A,+,\cdot_K)$ es un $K$-espacio vectorial
3. $(\lambda\cdot_K X)\cdot Y=\lambda\cdot_K(X\cdot Y)=X\cdot(\lambda\cdot_K Y)\quad\forall\lambda\in K,\forall X,Y\in A$
Se dice entonces que $A$ es una $K$-álgebra.
#### El producto de matrices describe un sistema lineal $n\times m$
El producto de matrices nos permite escribir un sistema lineal de $n$ ecuaciones con $m$ incógnitas
$$\begin{cases}a_{11}x_1+a_{12}x_2+\cdots+a_{1m}x_m & = &b_1\\  & \vdots & \\ a_{n1}x_1+a_{n2}x_2+\cdots+a_{nm}x_m & = & b_n\end{cases}$$
en la forma
$$A\cdot x=b$$
donde $A\in K^{n\times m}$ es la matriz asociada al sistema, $x\in K^{m\times 1}$ se define como $x_{i1}=x_i$, y $b\in K^{n\times 1}$ se define como $b_{j1}=b_j$. De esta manera, un sistema lineal puede verse como una única ecuación con una única incógnita $x$, pero que involucra matrices en lugares de escalares.
### Conceptos de matrices
#### Matrices triangulares
>[!note] Matriz triangular superior e inferior
>Una matriz $n\times n$ **triangular superior** $U=[u_{ij}]$ tiene, para cada $j=1,2,...,n$, las entradas $u_{ij}=0$, para cada $i=j+1,j+2,...,n$; y una matriz **triangular inferior** $L=[l_{ij}]$ tiene, para cada $j=1,2,...,n$ las entradas $l_{ij}=0$, para cada $i=1,2...,j-1$.

En consecuencia a esta definición, tenemos que una matriz diagonal es tanto triangular superior como triangular inferior debido a que sus entradas diferentes de cero deben estar en la diagonal principal.
#### Matriz transpuesta
>[!note] Matriz transpuesta
>Sea $A\in K^{n\times m}$. La matriz transpuesta de $A$, notada como $A^t$, es la matriz $A^t\in K^{m\times n}$ definida por $(A^t)_{ij}=A_{ji}$ para cada $1\leq i\leq m, 1\leq j\leq n$.
#### Traza de la matriz
>[!note] Traza de una matriz
>Sea $A\in K^{n\times n}$. Se llama **traza de la matriz** $A$, y se nota $\text{tr}(A)$, al escalar $\text{tr}(A)=\sum_{i=1}^n A_{ii}$.
### Matrices invertibles
No es cierto que todo elemento no nulo de $K^{n\times n}$ tenga inverso con respecto al producto.
>[!note] Matriz invertible
>Una matriz $A\in K^{n\times n}$ se dice invertible (no singular) si existe una matriz $B\in K^{n\times n}$ tal que $A\cdot B=B\cdot A=I_n$. 
#### Unicidad de la inversa de una matriz
La matriz $B$ es única. Puesto que si $A\cdot B=B\cdot A=I_n$ y $A\cdot C=C\cdot A=I_n$, entonces
$$B=I_n\cdot B=(C\cdot A)\cdot B=C\cdot(A\cdot B)=C\cdot I_n=C$$
#### Grupo lineal general
>[!note] Grupo lineal general
>$(\text{GL}(n,K),\cdot)$ es un grupo, llamado **grupo lineal general** $(n,K)$, el cual corresponde al conjunto de todas las matrices invertibles de tamaño $n\times n$ sobre un campo dado.
#### Propiedades de las matrices invertibles
Para cada $n\in\mathbb{N}$ consideraremos el conjunto de todas las matrices invertibles en $K^{n\times n}$:
$$\text{GL}(n,K)=\{A\in K^{n\times n}:A\ \text{es invertible}\}$$
Para cada $n\in\mathbb{N}$, se verifican las siguientes propiedades:
1. Si $A,B\in\text{GL}(n,K)$, entonces $A\cdot B\in\text{GL}(n,K)$. Más aún, $(A\cdot B)^{-1}=B^{-1}A^{-1}$. En particular, el producto de matrices $\cdot$ es una operación en $\text{GL}(n,K)$.
2. $I_n\in\text{GL}(n,K)$
3. Si $A\in\text{GL}(n,K)$, entonces $A^{-1}\in\text{GL}(n,K)$.
#### Estudio de invertibilidad de matrices
Para determinar si una matriz es invertible y hallar su inversa analizamos los siguientes puntos:
- ==Dada $A\in K^{n\times n}$, se arma una matriz ampliada en $K^{n\times 2n}$ cuyas primeras $n$ columnas corresponden a la matriz $A$ y cuyas últimas $n$ columnas están formadas por los $n$ vectores de la base canónica de $K^n$.== O sea, planteamos la ecuación $A\cdot B=I_n$ con $B\in K^{n\times n}$, la subdividimos en $n$ sistemas lineales $A\cdot B_i=e_i,1\leq i\leq n$, igualando columna a columna, y escribir la matriz ampliada de los $n$ sistemas.
- Si al triangular la matriz, no aparece **ningún cero en la diagonal**, se puede resolver los sistemas que resultan (que tienen solución única) y hallar entonces la inversa de $A$.
  ==Si al triangular la matriz aparece un cero en la diagonal, la matriz $A$ no es invertible.== Su aparición implica que el sistema homogéneo cuya matriz es $A$ tiene solución no trivial, es decir, existe $x_0\in K^n$ no nulo tal que $A\cdot x_0=0$. Si $A$ fuese invertible, multiplicando por $A^{-1}$ resultaría $x_0=A^{-1}\cdot A\cdot x_0=0$, contradiciendo que $x_0\neq0$. 
Estos puntos constituyen corolarios del [[Teórica 1 - Resolución de sistemas de ecuaciones lineales#Teorema de estructura para $ text{GL}(n,K)$|Teorema de estructura para el grupo lineal general]].
### Matrices elementales: Operaciones elementales como TL
>[!note] Matrices elementales
>Tipo especial de matrices que se obtienen a partir de la matriz identidad aplicando solo una operación elemental de fila o columna.

Existen **3 familias de matrices elementales**, correspondientes a los tres tipos de operaciones elementales permitidas (Gauss).
#### Intercambio de filas o ecuaciones como TL
Sean $1\leq i,j\leq n$. Se define $P^{ij}\in K^{n\times n}$ como
$$P^{ij}=I_n-E^{ii}-E^{jj}+E^{ij}+E^{ji}$$
Observamos que $P^{ij}$ es la matriz que resulta al intercambiar las filas $i$ y $j$ en la matriz $I_n$.
Es fácil verificar que, dada $B\in K^{n\times n}$ el producto $P^{ij}B$ es la matriz que resulta al intercambiar en la matriz $B$ las filas $i$ y $j$.
En particular, $P^{ij}\cdot P^{ij}=I_n$, es decir que $P^{ij}$ es invertible y $P^{ij^{-1}}=P^{ij}$.
#### Multiplicación fila por escalar no nulo
Sea $a\in K,a\neq 0$, y sea $1\leq i\leq n$. Se define $M_i(a)\in K^{n\times n}$ como
$$M_i(a)=I_n+(a-1)\cdot E^{ii}$$
Siendo $M_i(a)$ la matriz que se obtiene al multiplicar por $a$ la $i$-ésima fila de la matriz $I_n$.
Dada $B\in K^{n\times n}$ se tiene que
$$(M_i(a)\cdot B)_{kj}=\begin{cases}B_{kj} & \text{si }k\neq j  \\ a\cdot B_{kj} & \text{si }k=j\end{cases}$$
O sea, $M_i(a)\cdot B$ es la matriz que resulta al multiplicar por $a$ la $i$-ésima fila de $B$.
En particular, $M_i(a)\cdot M_i(a^{-1})=M_i(a^{-1})\cdot M_i(a)=I_n$, de donde $M_i(a)\in\text{GL}(n,K)$ y $(M_i(a))^{-1}=M_i(a^{-1})$.
#### Reemplazar una ecuación por combinación lineal con otras
Sea $a\in K$ y sean $i\neq j$ con $1\neq i,j\neq n$. Se define $T^{ij}(a)\in K^{n\times n}$ como
$$T^{ij}(a)=I_n+a\cdot E^{ij}$$
La matriz que se obtiene de la matriz $I_n$ al sumarle a la $i$-ésima fila, $a$ por la fila $j$.
Si $B\in K^{n\times n}$, entonces
$$(T^{ij}(a)\cdot B)_{kl}=\begin{cases}B_{kl} & \text{si }k\neq i \\ B_{il}+a\cdot B_{jl} & \text{si }k=i\end{cases}$$
O sea que $T^{ij}(a)\cdot B$ es la matriz que se obtiene de $B$ al sumarle a la $i$-ésima fila, $a$ por la fila $j$.
En particular, se tiene que $T^{ij}(a)\cdot T^{ij}(-a)=T^{ij}(-a)\cdot T^{ij}(a)=I_n$, con lo que $T^{ij}(a)\in\text{GL}(n,K)$ y $(T^{ij}(a))^{-1}=T^{ij}(-a)$.
### Observación teórica de la obtención de una matriz inversa
Sea $A\in K^{n\times n}$. Entonces existen matrices elementales $E_1,...,E_r\in K^{n\times n}$ tales que $E_r ... E_1\cdot A$ es triangular superior. Si además, $E_r....E_1\cdot A$ no tiene ceros en la diagonal, existen matrices elementales $E_{r+1},...,E_s$ tales que $E_s...E_{r+1}\cdot E_r...E_1\cdot A=I_n$. En consecuencia, $A$ es producto de matrices elementales: $A=E_1^{-1}...E_s^{-1}$, y $A^{-1}=E_s...E_1$.
O sea, si mediante la aplicación de operaciones elementales a las filas de la matriz $A$ obtenemos la matriz identidad $I$, entonces aplicando las mismas operaciones en las filas de $I$ obtendremos $A^{-1}$.
### Teorema de estructura para $\text{GL}(n,K)$
>[!note] Teorema de estructura para $\text{GL}(n,K)$
>Así como el Teorema Fundamental de la Aritmética en $\mathbb{Z}$ dice que todo número entero no nulo es producto de enteros primos, la observación anterior nos dice que ==toda matriz en $\text{GL}(n,K)$ es producto de matrices elementales.==
### Coordenadas

#### Coordenadas de un vector en una base
Dado un $K$-espacio vectorial $V$ de dimensión $n$ y fijada una base $B$ en $V$, mediante el concepto de coordenadas de un vector en la base $B$ podremos identificar cada elemento de $V$ con un vector $K^n$ y trabajar entonces con elementos de $K^n$.
>[!note] Vector de coordenadas de $x$ en la base $B$
>Sea $V$ un $K$-espacio vectorial de dimensión finita y sea $B=\{v_1,...,v_n\}$ una base de $V$. Dado $x\in V$, existen únicos $\alpha_1,...,\alpha_n\in K$ tales que $x=\alpha_1 v_1+...+\alpha_n v_n$. El vector $(\alpha_1,...,\alpha_n)\in K^n$ es el **vector de coordenadas** de $x$ en la base $B$ y será denotado por $(x)_B$.
#### Cambios de base
Dadas dos bases de un mismo $K$-espacio vectorial $V$ de dimensión finita, cada elemento de $V$ tiene asociados dos vectores de coordenadas, uno en cada una de las bases. Mediante la **matriz de cambio de base** podemos obtener las coordenadas de un vector con respecto a una base de $V$ a partir de las coordenadas del vector en otra base.
>[!note] Matriz de cambio de base
>Sea $V$ un $K$-espacio vectorial de dimensión $n$, y sean $B_1=\{v_1,...,v_n\}$ y $B_2=\{w_1,...,w_n\}$ bases de $V$. Para cada $1\leq j\leq n$, sean $\alpha_{ij}\in K\quad(1\leq i\leq n)$ tales que $v_j=\sum_{i=1}^n\alpha_{ij}w_i$.
>Definimos a la matriz de cambio de base de $B_1$ a $B_2$, y se nota $C(B_1,B_2)\in K^{n\times n}$, a la matriz definida por $(C(B_1,B_2))_{ij}=\alpha_{ij}$ para cada $1\leq i,j\leq n$.
##### Propiedad
Sea $V$ un $K$-espacio vectorial de dimensión finita, y sean $B_1$ y $B_2$  bases de $V$. Entonces, para cada $x\in V$,
$$C(B_1,B_2)\cdot((x)_{B_1})^t=((x)_{B_2})^t$$
### Determinante de una matriz
El **determinante** de una matriz ==proporciona resultados de existencia y unicidad== para sistemas lineales que tienen el mismo número de incógnitas y de ecuaciones. Es denotado como $\text{det}A$ o como $|A|$.
>[!note] Determinante de una matriz
>Sea $A\in K^{n\times n}$.
>1. Si $A=[a]$ es una matriz $1\times 1$, entonces $\text{det}A=a$.
>2. Si $A$ es una matriz $n\times n$, con $n>1$, el menor $M_{ij}$ es el determinante de la submatriz $(n-1)\times(n-1)$ de $A$ obtenida al quitar la $i$-ésima fila y la $j$-ésima columna de la matriz $A$.
>3. El **cofactor** $A_{ij}$ asociado con $M_{ij}$ está definido por $A_{ij}=(-1)^{i+j}M_{ij}$.
>4. El **determinante** de la matriz $A$, cuando $n>1$, está dado ya sea por $$\text{det}A=\sum_{j=1}^na_{ij}A_{ij}=\sum_{j=1}^n(-1)^{i+j}a_{ij}M_{ij},\qquad\text{para }i=1,2,...,n$$o mediante $$\text{det}A=\sum_{i=1}^n(-1)^{i+j}a_{ij}M_{ij},\qquad\text{para cualquier }j=1,2,...,n$$
#### Complejidad algorítmica
Notar que esta definición es recursiva, de modo que calcular el determinante de una matriz general $n\times n$ requiere $O(n!)$ multiplicaciones/divisiones y sumas/restas. Incluso para valores relativamente pequeños de $n$, el número de cálculos se vuelve difícil de manejar.
>[!tip] Solución óptima: fila/columna con mayor cantidad de ceros
>Sin embargo, como existen $2n$ definiciones diferentes de $\text{det}A$, dependiendo de la fila o columna seleccionada, todas las definiciones llevan al mismo resultado numérico. La flexibilidad en la definición se utiliza en le siguiente ejemplo. Es más conveniente calcular $\text{det}A$ a lo largo de la fila o de la columna con la mayor cantidad de ceros.
#### Propiedades
Suponga que $A$ es una matriz $n \times n$:
1. Si cualquier fila o columna $A$ solo tiene entradas cero, entonces $\det A = 0$.
2. Si $A$ tiene dos filas o dos columnas iguales, entonces $\det A = 0$.
3. Si $\tilde{A}$ se obtiene a partir de $A$ mediante la operación $(E_i) \leftrightarrow (E_j)$, con $i \neq j$, entonces $\det \tilde{A} = -\det A$.
4. Si $\tilde{A}$ se obtiene a partir de $A$ mediante la operación $(\lambda E_i) \rightarrow (E_i)$, entonces $\tilde{A} = \lambda \det A$.
5. Si $\tilde{A}$ se obtiene a partir de $A$ mediante la operación $(E_i + \lambda E_j) \rightarrow (E_i)$ con $i \neq j$, entonces $\det \tilde{A} = \det A$.
6. Si $B$ también es una matriz $n \times n$, entonces $\det AB = \det A \det B$.
7. $\det A^t = \det A$.
8. Cuando $A^{-1}$ existe, $\det A^{-1} = (\det A)^{-1}$.
9. Si $A$ es una matriz triangular superior, triangular inferior o diagonal, entonces $A = \prod_{i=1}^{n} a_{ii}$.
Resumimos estas propiedades en la siguiente tabla
![[Pasted image 20240812235536.png]]
#### Declaraciones equivalentes para matrices cuadradas
Las siguientes declaraciones son equivalentes para cualquier matriz $A\in K^{n\times n}$.
1. La ecuación $A\,x=0$ tiene solución única $x=0$.
2. El sistema $A\,x=b$ tiene solución única para cualquier vector de columna $n$-dimensional $b$.
3. La matriz $A$ no es singular, es decir, existen $A^{-1}$.
4. $\text{det}A\neq0$
5. La eliminación gaussiana con intercambios de fila se puede realizar en el sistema $A\,x=b$ para cualquier vector de columna $n$-dimensional $b$.
# Sistemas de ecuaciones lineales
## Sistemas lineales homogéneos
>[!note] Sistema lineal homogéneo
>Un sistema lineal homogéneo de $n$ ecuaciones con $m$ incógnitas a coeficientes en un cuerpo $K$ es un sistema del tipo:
>$$\begin{cases}a_{11}x_1+a_{12}x_2+\cdots+a_{1m}x_m=0 \\ a_{n1}x_1+a_{n2}x_2+\cdots+a_{nm}x_m=0\end{cases}$$
>donde $a_{ij}\in K$ para cada $1\leq i\leq n,1\leq j\leq m$.
### Soluciones de un sistema lineal homogéneo
>[!tip] Conjunto de soluciones de un sistema lineal homogéneo
>El conjunto de las soluciones de un sistema lineal homogéneo con $m$ incógnitas es un subespacio de $K^m$. O sea, resolver este tipo de sistemas nos da el **sistema de generadores** para el subespacio de las soluciones.

Para resolver estos sistemas debemos realizar una serie de operaciones que nos permitan obtener otro sistema, pero con el mismo conjunto de soluciones. Estos sistemas son llamados **sistemas equivalentes**.
### Sistemas equivalentes
>[!note] Sistemas equivalentes
>Dos sistemas lineales homogéneos se dicen **equivalentes** si sus conjuntos de soluciones son iguales.

Por ejemplo: los siguientes sistemas lineales homogéneos a coeficientes en $\mathbb{R}$ son equivalentes:
$$\begin{cases}x+y+z&=0\\y+z&=0\end{cases}\qquad\begin{cases}x&=0\\y+z&=0\end{cases}$$
## Triangulación: eliminación gaussiana
>[!note] Eliminación gaussiana
>Método eficiente para la resolución de sistemas de ecuaciones lineales en forma directa (es decir, sin utilizar métodos iterativos que aproximan la solución).

Un ejemplo sería la resolución del siguiente sistema de ecuaciones.
$\begin{cases} x+5y+5z = 2 \\ 2x+2y-3z = -1 \\ -x-9y+2z = 9 \end{cases}$
A partir de este sistema de ecuaciones, construimos la **matriz ampliada** de coeficientes y términos independientes:
$$\left
(\begin{array}{ccc|c}
1 & 5 & 5 & 2 \\
2 & 2 & -3 & -1 \\
-1 & -9 & 2 & 9
\end{array}\right)
$$
Triangulamos la matriz realizando operaciones de filas. Las operaciones permitidas (que no afectan las soluciones del sistema) son:
- Sumar (o restar) a una fila un múltiplo de otra fila,
- intercambiar dos filas entre sí,
- multiplicar una fila por un escalar distinto de 0.
Entonces la resolución de este sistema sería el siguiente:
$$\begin{aligned}
	&\begin{pmatrix}
		1 & 5 & 5 &|& 2\\
		2 & 2 & -3 &|& -1\\
		-1 & -9 & 2 &|& 9
	\end{pmatrix}\xrightarrow{f_2-2f_1\to f_2}\begin{pmatrix}
		1 & 5 & 5 &|& 2\\
		0 & -8 & -13 &|& -5\\
		-1 & -9 & 2 &|& 9
	\end{pmatrix}\xrightarrow{f_3+f_1\to f_3}\\
	&\begin{pmatrix}
		1 & 5 & 5 &|& 2\\
		0 & -8 & -13 &|& -5\\
		0 & -4 & 7 &|& 11
	\end{pmatrix}\xrightarrow{f_3-\frac12f_2\to f_3}\begin{pmatrix}
		1 & 5 & 5 &|& 2\\
		0 & -8 & -13 &|& -5\\
		0 & 0 & \frac{27}2 &|& \frac{27}2
	\end{pmatrix}
\end{aligned}$$

#### Sustitución hacia atrás (Back-substitution):
Del sistema triangular superior equivalente:
$$\begin{aligned}
\frac{27}{2}z = \frac{27}{2} &\implies z = 1 \\
-8y - 13(1) = -5 &\implies -8y = 8 \implies y = -1 \\
x + 5(-1) + 5(1) = 2 &\implies x = 2
\end{aligned}$$
Solución única: $(x, y, z) = (2, -1, 1)$.

>[!tip] EliminacionGaussianaYSustitucionAtras
>**Cuándo usarla**: Para resolver sistemas lineales cuadrados o rectangulares $Ax=b$ de forma analítica o directa sin invertir $A$.
>**Pasos**:
>1. Construir la matriz ampliada $(A \mid b)$.
>2. Para cada columna $k = 1, \dots, n-1$, elegir pivote no nulo $a_{kk} \neq 0$ (si es cero, permutar filas $P_{ik}$).
>3. Calcular multiplicadores $m_{ik} = \frac{a_{ik}}{a_{kk}}$ y restar $f_i \leftarrow f_i - m_{ik} f_k$ para todo $i > k$ hasta obtener la matriz triangular superior $(U \mid y)$.
>4. Aplicar sustitución hacia atrás desde la última incógnita hacia arriba: $x_i = \frac{1}{u_{ii}}\left(y_i - \sum_{j=i+1}^n u_{ij}x_j\right)$.
>**Error típico**: Confundir el orden del multiplicador (dividir por el elemento a anular en lugar de dividir por el pivote) o aplicar la operación sólo sobre la matriz $A$ olvidando actualizar el vector término independiente $b$.

**Referencias y Conexiones**:
- *Se usa en*: [[Práctica 1 - Resolución de sistemas de ecuaciones lineales]] Ej. 1 y Ej. 2.
- *Fundamento de*: [[Teórica 3 - Sistemas lineales y factorización#Factorización LU]], donde los multiplicadores $m_{ik}$ forman la matriz triangular inferior $L$.
### Cantidad de soluciones de un sistema homogéneo
#### Teorema de existencia de soluciones no triviales en sistemas homogéneos donde $n<m$
>[!note] Soluciones no triviales de sistemas homogéneos donde $n<m$
>Sea $H$ un sistema lineal homogéneo de $n$ ecuaciones con $m$ incógnitas. Supongamos que $n<m$. Entonces existe $x\in K^m,x\neq0$, que es solución del sistema $H$.
##### Demostración
Por **inducción** en la cantidad $n$ de ecuaciones de $H$.
Si $n=1,m\geq2$, entonces $H:a_{11}x_1+a_{12}x_2+\cdots+a_{1m}x_m=0$.
- Si $a_{11}=0$, entonces (1,0,...,0) es solución del sistema
- Si $a_{11}\neq0$, entonces $(\frac{-a_{12}}{a_{11}},1,0,...,0)$ es solución.
Supongamos que el resultado vale para sistemas con $n$ ecuaciones y sea $H$ un sistema de $n+1$ ecuaciones con $m$ incógnitas, $n+1<m$.
Triangulando la matriz del sistema, resulta que es equivalente a una de la forma:
$$\begin{pmatrix}a_{11} & a_{12} & \cdots & a_{1m}\\0 & & B\end{pmatrix}$$
donde $B\in K^{n\times(m-1)}$, y $m-1>n$.
Por lo tanto, el sistema cuya matriz asociada es $B$ está en las condiciones de la hipótesis inductiva. Luego, existe $(x_1,...,x_{m-1})\neq0$ que es la solución del sistema asociado a $B$.
- Si $a_{11}=0$, entonces $(1,0,...,0)$ es solución del sistema original.
- Si $a_{11}\neq0$, entonces $(-\frac{1}{a_{11}}\cdot(\sum_{i=2}^ma_{1i}x_{i-1}),x_1,...,x_{m-1})$ es una solución no nula del sistema.
#### Teorema de existencia de soluciones no triviales en sistemas homogéneos donde $n=m$
>[!note] Soluciones no triviales de sistemas homogéneos de $n$ incógnitas y ecuaciones
>Sea $H$ un sistema lineal homogéneo de $n$ ecuaciones y $n$ incógnitas. Sea $H'$ el sistema equivalente a $H$ cuya matriz $B$ es triangular superior. Entonces,
>$$H \text{es SCD} \iff B_{ii}\neq0,\forall 1\leq i\leq n$$
##### Demostración
###### $\Longleftarrow$
Supongamos que $B=\begin{pmatrix}B_{11} & \cdots & B_{1n} \\ 0 & \ddots & \vdots \\ \cdots & 0 & B_{nn}\end{pmatrix}$ con $B_{ii}\neq0,\forall 1\leq i\leq n$.
Entonces, la última ecuación del sistema $H'$ es $B_{nn}x_n=0$ y, como $B_{nn}\neq0$, resulta que $x_n=0$. Reemplazando en la ecuación $x_n$ por 0, queda $B_{n-1\ n-1}x_{n-1}=0$, de donde $x_{n-1}=0$.
Siguiendo de este modo, para cada $k=n-2,...,1$ de la $k$-ésima ecuación se obtiene $x_k=0$.
###### $\implies$
Supongamos que el sistema homogéneo $H$ (equivalente a $H'$) es **SCD**, es decir, tiene **única solución** (la trivial).
Sea $B$ su matriz triangular superior. Supongamos, por el contrario, que existe algún índice $i$ tal que: $B_{11}\neq0,...,B_{ii}\neq0$ y $B_{i+1\ i+1}=0$, o sea
$$\begin{pmatrix}
	B_{11} & \cdots & B_{1\,i} & B_{1\,i+1} & \cdots & B_{1n}\\
	0 & \ddots & \vdots & \vdots & & \vdots\\
	0 & \cdots & B_{ii} & B_{i\,i+1} & \cdots & B_{i\,n}\\
	0 & \cdots & 0 & 0 & \cdots & B_{i+1,n}\\
	\vdots & & \vdots & \vdots & \ddots & \vdots\\
	0 & \cdots & 0 & 0 & \cdots & B_{nn}
\end{pmatrix}$$
Entonces, el sistema $H'$ tiene la forma:
$$\begin{cases}
	B_{11}x_1+\cdots+B_{1n}x_n=0\\
	\vdots\\
	0\cdot B_{ii}x_i+\cdots+B_{in}x_n=0\\
	\vdots\\
	0\cdot x_n=0
\end{cases}$$
Observá que en la ecuación $i+1$-iésima **no aparece** $x_{i+1}$ (porque $B_{i+1\,i+1}=0$).
Puesto que la variable $x_{i+1}$ no está determinada por el sistema, queda como un parámetro libre. De este modo, admite soluciones no triviales.
Por lo tanto, por reducción al absurdo, concluimos que si el sistema es SCD, entonces $B_{ii}\neq0,\forall i\in[1,n]$.
## Matriz escalonada
>[!note] Matriz escalonada
>Una matriz $A\in \mathbb{K}^{m\times n}$ es escalonada si para todo $2\leq i \leq m$, la fila $i$ tiene más ceros iniciales que la fila $i-1$.

Por ejemplo:
- $A=\begin{pmatrix} 1 & 2 & 0 & 7 & -1 \\ 0 & 1 & 2 & 4 & -1 \\ 0 & 0 & 0 & 0 & 5\end{pmatrix}$ es una matriz escalonada
- $A=\begin{pmatrix} 1 & 2 & 0 & 7 & -1 \\ 0 & 0 & 0 & 4 & -1 \\ 0 & 0 & 0 & -3 & 5\end{pmatrix}$ **no** es una matriz escalonada
Si a partir de la matriz ampliada, obtenemos una matriz escalonada equivalente, podemos resolver el sistema por sustitución hacia atrás, despejando en cada fila la primera variable con coeficiente no nulo en función de las siguientes variables.
### Algoritmo de eliminación gaussiana
Si consideramos que el *input* es una matriz ampliada $A\in\mathbb{R}^{m\times(n+1)}$, y el output es una matriz escalonada equivalente.
Para $i=1,2,...,m-1$:
1. Consideramos la matriz $\tilde{A}$ tomando las filas de $A$ desde la fila $i$ hasta la última.
2. Buscamos $k$ tal que la columna $k$ de $\tilde{A}$ tenga elementos no nulos, y todas las columnas anteriores sean nulas.
3. Para $j=i+1,...,m$, calculamos el multiplicador $m_{ji} = \frac{a_{jk}}{a_{ik}}$ y realizamos en $A$ la operación elemental:
$$f_j - m_{ji} f_i \to f_j$$

#### Análisis de Costo Computacional (Flops)
Para una matriz cuadrada de $n \times n$:
- En el paso $k$ ($k=1,\dots,n-1$), se calculan $n-k$ multiplicadores (divisiones).
- Para cada una de las $n-k$ filas inferiores, se actualizan sus $n-k$ columnas no nulas, realizando 1 multiplicación y 1 resta por entrada: $2(n-k)^2$ flops.
Sumando sobre todas las etapas de eliminación:
$$\text{Flops de Triangulación} = \sum_{k=1}^{n-1} 2(n-k)^2 \approx 2 \int_0^n (n-x)^2 dx = \frac{2}{3}n^3 \text{ flops}$$
- La actualización del lado derecho $b$ toma $\sum_{k=1}^{n-1} 2(n-k) \approx n^2$ flops.
- La **sustitución hacia atrás** requiere: $\frac{n(n+1)}{2}$ productos/restas y $n$ divisiones $\approx n^2$ flops.

#### Consecuencias
1. **Dominancia asintótica**: La eliminación gaussiana está dominada en un 99% por la triangulación ($\frac{2}{3}n^3$ flops). Resolver el sistema triangular posterior es computacionalmente despreciable ($O(n^2)$).
2. **Impacto en Factorización**: Si se necesita resolver $Ax = b$ para múltiples vectores $b$, recalcular Gauss costaría $k \times \frac{2}{3}n^3$. Esto motiva la factorización $LU$ de la [[Teórica 3 - Sistemas lineales y factorización#Factorización LU]], que factoriza $A$ una sola vez ($\frac{2}{3}n^3$) y luego resuelve cada $b$ en solo $2n^2$ flops.
### Matriz escalonada en Python
Como en los paquetes comunes de **Python** no existe un comando para eliminación gaussiana, utilizamos el siguiente programa que desarrolla una función llamada `row_echelon` (nombre en inglés para una matriz escalonada por filas).
``` Python
import numpy as np

def matriz_escalonada(A):
	""" Devuelve la Forma Escalonada de la matriz A """
	A = np.copy(M)
	if (issubclass(A.dtype.type, np.integer)):
		A = A.astype(float)
		
	# si la matriz A no tiene columnas o filas,
	# ya está en REF (Row Echelon Form), por lo que devolvemos la misma
	r, c = A.shape
	if r == 0 or c == 0:
		return A
		
	# buscamos un elemento no nulo en la primera columna
	for i in range(len(A)):
		if A[i,0] != 0:
			break
			
	else:
		# si todos los elementos en la primera columna son cero,
		# realizamos REF en la matriz a partir de la segunda columna
		B = row_echelon(A[:,1:])
		# y luego añadimos la primera columna de ceros de vuelta
		return np.hstack([A[:,:1], B])
		
	# si el elemento no nulo no está en la primera fila,
	# intercambiamos las filas
	if i > 0:
		ith_row = A[i].copy()
		A[i] = A[0]
		A[0] = ith_row
	
	# dividimos la primera fila por el primer elemento en ella
	A[0] = A[0] / A[0,0]
	# restamos todas las filas subsiguientes con la primera fila 
	# (ahora tiene 1 como primer elemento) multiplicado por el
	# elemento correspondiente en la primera columna
	A[1:] -= A[0] * A[1:,0:1]
	
	# realizamos REF en la matriz a partir de la segunda fila, desde la
	# segunda columna
	B = row_echelon(A[1:,1:])
	
	# añadimos la primera fila y la primera columna (cero), y devolvemos
	return np.vstack([A[:1], np.hstack([A[1:,:1],B])])
```
## Clasificación de sistemas de ecuaciones
Los sistemas de ecuaciones pueden ser clasificadas según cuántas soluciones tiene. Estas soluciones pueden ser determinadas a partir de la matriz ampliada.
Si al escalonar la matriz ampliada, llegamos a una ecuación:
$$0x_1+\cdots+0x_n=a\neq0$$
el sistema no tiene solución. Por el contrario, si en todas las filas que se anulan los coeficientes de las variables, obtenemos también 0 en el término independiente, el sistema admite solución (podemos ir resolviendo hacia atrás).
De esta forma, obtenemos los siguientes casos para un sistema de $m$ ecuaciones y $n$ incógnitas.
- **Sistema Incompatible:** Sistema que no tiene solución. Al escalonar obtenemos una ecuación $0x_1+\cdots+0x_n=a\neq0$. Por ejemplo: $$\left(\begin{array}{ccc|c}5 & 3 & 2 & 11 \\ 0 & 9 & -1 & 10 \\ 0 & 0 & 0 & 4 \end{array}\right)$$
- **Sistema Compatible Determinado:** Sistema con solución única. Al escalonar la matriz ampliada, obtenemos exactamente $n$ ecuaciones no nulas, y podemos obtener la solución despejando las variables. Por ejemplo: $$\left(\begin{array}{ccc|c}5 & 3 & 2 & 11 \\ 0 & 9 & -1 & 10 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{array}\right)$$
- **Sistema Compatible Indeterminado:** Sistema con infinitas soluciones. Al escalonar la matriz ampliada obtenemos menos de $n$ filas no nulas en las primeras $n$ columnas, y el resto de las filas son nulas en todas las columnas. Por ejemplo: $$\left(\begin{array}{ccc|c}5 & 3 & 2 & 11 \\ 0 & 9 & -1 & 10 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{array}\right)$$ 
# Espacios vectoriales y bases
>[!note] Espacio vectorial
>Conjunto de vectores que constituyen la **estructura base de los conjuntos**, formados de tal manera que *es posible realizar la suma y la multiplicación* entre sus elementos (podemos generar combinaciones lineales de estos).
## Espacios vectoriales
Un espacio vectorial requiere la existencia de dos conjuntos:
- Conjunto $K$ (escalares)
- Conjunto $V$ (vectores)
Ambos conjuntos deben satisfacer ciertas propiedades, que esencialmente se refieren a que los elementos de $V$ se puedan sumar entre sí y multiplicar por elementos de $K$.
### Definiciones base para un espacio vectorial
#### Operación (binaria o ley de composición interna)
>[!note] Operación binaria
>Sea $A$ un conjunto no vacío. Definimos a la *operación* de $A$ como una función $*:A\times A\to A$.
>*Notación*: $*(a,b)=c$ se escribe $a*b=c$.
##### Ejemplos
- $+:\mathbb{N}\times\mathbb{N}\to\mathbb{N}$, tal que $+(a,b)=a+b$, es una operación de $\mathbb{N}$
- Como la resta, $-(a,b)=a-b$, no es una función de $\mathbb{N}\times\mathbb{N}\to\mathbb{N}$, entonces no es una operación $\mathbb{N}$
- La suma $+$, el producto $\cdot$ y la resta $-$ son operaciones de $\mathbb{Z},\mathbb{Q},\mathbb{R},\mathbb{C}$
##### Propiedades básicas de la operación binaria
1. $*$ se dice **asociativa** si $(a*b)*c=a*(b*c)\ \forall a,b,c\in A$.
2. Se dice que $*$ tiene *elemento neutro* si $\exists e\in A$ tal que $e*a=a*e=a$ para cada $a\in A$.
   %% Notar que si $*$ tiene elemento neutro, éste es ÚNICO, ya que si $e$ y $e'$ son elementos neutros, $e'=e*e'=e$ %%
3. Si $*$ tiene elemento neutro $e$, se dice que todo elemento tiene *inverso* para $*$ si $\forall a\in A, \exists a'\in A$ tal que $a*a'=a'*a=e$.
4. Se dice **conmutativa** si $a*b=b*a\ \forall a,b\in A$.
#### Grupo y grupo abeliano (conmutativo)
>[!note] Grupo y grupo abeliano
>Sea $A$ un conjunto, y sea $*$ una operación en $A$ que satisface las propiedades **1.**, **2.** y **3.** de la definición anterior. Entonces $(A,*)$ se llama *grupo*. Si además $*$ cumple la propiedad **4.**, se dice que $(A,*)$ es un *grupo abeliano* o *conmutativo*.
##### Ejemplos
- $(\mathbb{N},+)$ no es un grupo: se puede probar que no tiene elemento neutro.
- $(\mathbb{Z},+),(\mathbb{Q},+),(\mathbb{R},+)$ y $(\mathbb{C},+)$ son grupos abelianos.
- $(\mathbb{Z},\cdot)$ no es un grupo: se puede probar que sólo $1$ y $-1$ tienen inverso multiplicativo.
- $(\mathbb{Q}-{0},\cdot),(\mathbb{R}-{0},\cdot),(\mathbb{C}-{0},\cdot)$ son grupos abelianos
- $A=\{f:\mathbb{R}\to\mathbb{R}\}, *=\circ$ (composición de funciones). Entonces $(A,*)$ no es un grupo: las únicas funciones con inversa para $\circ$ son las biyectivas.
- $S_{\mathbb{R}}=\{f:\mathbb{R}\to\mathbb{R}/f\ \text{es biyectiva}\}$, $*=\circ$. Entonces ($S_{\mathbb{R}},\circ)$ es un grupo.
- $C$ un conjunto, $\mathbb{P}(C)=\{S\subset C\}$. Se define como la operación $\Delta:\mathbb{P}(C)\times\mathbb{P}(C)\to\mathbb{P}(C)$, llamada diferencia simétrica de la siguiente forma: $$A\Delta B=(A\cup B)-(a\cap B)$$Entonces $(\mathbb{P}(C),\Delta)$ es un grupo abeliano.
#### Anillo
>[!note] Anillo
>Sea $A$ un conjunto y sean $+$ y $\cdot$ operaciones de $A$. Se dice que $(A,+,\cdot)$ es un *anillo* si:
>1. $(A,+)$ es un grupo abeliano.
>2. $\cdot$ es asociativa y tiene elemento neutro
>3. Valen las propiedades distributivas: Para $a,b,c\in A$,
>	- $a\cdot(b+c)=a\cdot b+ a\cdot c$
>	- $(b+c)\cdot a=b\cdot a+c\cdot a$
>Además, si $\cdot$ es conmutativa, se dice que $(A,+,\cdot)$ es un *anillo conmutativo*.
##### Ejemplos
- $(\mathbb{Z},+,\cdot),(\mathbb{Q},+,\cdot),(\mathbb{R},+,\cdot),(\mathbb{C},+,\cdot)$ son anillos conmutativos.
- $(\mathbb{Z}_n,+,\cdot)$ es un anillo conmutativo.
- Si $(A,+,\cdot)$ es un anillo conmutativo, entonces $(A[X],+,\cdot)$ es un anillo conmutativo con las operaciones usuales de polinomios.
- Si $C$ es un conjunto $(\mathbb{P}(C),\Delta,\cap)$ es un anillo conmutativo.
- $\{f:\mathbb{R}\to\mathbb{R}\}$ con las operaciones usuales de suma y producto de funciones es un anillo conmutativo.
#### Dominio de integridad
>[!note] Dominio de integridad
>Un anillo conmutativo $(A,+,\cdot)$ se llama un *dominio de integridad* o *dominio íntegro* si $a\cdot b=0 \implies a=0 \lor b=0$.
##### Ejemplos
- $(\mathbb{Z},+,\cdot),(\mathbb{Q},+,\cdot),(\mathbb{R},+,\cdot),(\mathbb{C},+,\cdot)$ son dominios de integridad.
- Si $A$ es un dominio de integridad, entonces $A[X]$ es un dominio de integridad.
- $\mathbb{Z}_p$ es un dominio de integridad $\iff$ $p$ es primo
#### Cuerpo
>[!note] Cuerpo
>Sea $K$ un conjunto, y sean $+$ y $\cdot$ operaciones de $K$. Se dice que $(K,+,\cdot)$ es un *cuerpo* si $(K,+,\cdot)$ es un anillo conmutativo y todo elemento no nulo de $K$ tiene inverso multiplicativo. Es decir:
>1. $(K,+)$ es un **grupo abeliano**
>2. $(K-{0},\cdot)$ es un **grupo abeliano**
>3. Vale la propiedad distributiva de $\cdot$ con respecto a $+$
##### Ejemplos
- $(\mathbb{Q},+,\cdot),(\mathbb{R},+,\cdot),(\mathbb{C},+,\cdot)$ son cuerpos
- $(\mathbb{Z}_p,+,\cdot)$ es un cuerpo $\iff$ $p$ es primo
#### Acción
Clase especial de funciones que se aplican a elementos de dos conjuntos distintos.
>[!note] Acción
>Sean $A$ y $B$ dos conjuntos. Una *acción* de $A$ en $B$ es una función $\cdot:A\times B\to B$.
>Notación: $\cdot(a,b)=a\cdot b$
### Espacio vectorial
>[!note] Espacio vectorial
>Sea $(K,+,\cdot)$ un cuerpo. Sea $V$ un conjunto no vacío, sea $+$ una operación en $V$ y sea $\cdot$ una acción de $K$ en $V$. Se dice que $(V,+,\cdot)$ es un *$K$-espacio vectorial* si se cumplen las siguientes condiciones:
>1. $(V,+)$ es un **grupo abeliano** (cerrado para la suma)
>2. La acción $\cdot:K\times V\to V$ satisface: (producto por escalar)
>	- $a\cdot(v+w)=a\cdot v+a\cdot w\quad \forall a\in K; \forall v,w\in V$
>	- $(a+b)\cdot v= a\cdot v+b\cdot v\quad \forall a,b\in K;\forall v\in V$
>	- $1\cdot v=v\quad v\in V$
>	- $(a\cdot b)\cdot v= a\cdot(b\cdot v)\quad \forall a,b\in K;\forall v\in V$

Los elementos de $V$ se llaman *vectores* y los elementos de $K$ se llaman *escalares*. La acción $\cdot$ se llama *producto por escalares*.
## Subespacios
Dentro de un $K$-espacio vectorial $V$, hay subconjuntos que heredan la estructura de $V$.
>[!note] Subespacios
>Sea $V$ un $K$-espacio vectorial. Un subconjunto $S\subset V$ **no vacío** se dice un *subespacio de* $V$ si la suma y el producto por escalares (de $V$) son una operación y una acción en $S$ que lo convierten en un $K$-espacio vectorial.

>[!quote] Definición de Strang
>Un subespacio vectorial es un conjunto no vacío que satisface los requisitos de un espacio vectorial, las combinaciones lineales permanecen en el subespacio. De este modo:
>- Si nosotros sumamos dos vectores $x$ e $y$ en el subespacio, $x+y$ pertenece al subespacio.
>- Si multiplicamos cualquier vector $x$ en el subespacio por un escalar $c$, $cx$ está en el subespacio.
### Verificación de un subespacio
Entonces, si queremos averiguar si un conjunto $S$ es un subespacio de $V$, éste lo será si y sólo si valen las siguientes condiciones:
1. $0\in S$. Verifica que $S$ no es un conjunto vacío. Entonces la podríamos reemplazar por $S\neq\emptyset$
2. **Cerrados para la suma:** $v,w\in S\implies v+w\in S$. Comprueba que $+$ es una operación de $S$
3. **Cerrados para el producto por escalares:** $\lambda\in K,v\in S\implies\lambda\cdot v\in S$. Se prueba que $\cdot$ es una acción
### Intersección de subespacios
>[!note] Intersección de subespacios
>Si $V$ es un $K$-espacio vectorial, $S$ y $T$ subespacios de $V$, entonces $S\cup T$ no es necesariamente un subespacio de $V$.
#### Ejemplo
Consideremos en $\mathbb{R}^2$ los subespacios $S=<(1,0)>$ y $T=<(0,1)>$. Observamos que $(1,0)\in S$ y $(0,1)\in T$; luego, ambos pertenecen a $S\cup T$. Pero, $(1,0)+(0,1)=(1,1)\not\in S\cup T$, puesto que $(1,1)\not\in S$ y $(1,1)\not\in T$.
### Espacio columna
>[!note] Espacio columna
>Espacio que contiene todas las combinaciones lineales de las columnas de $A$. Este corresponde a un subespacio de $K^m$.
#### Ejemplo
Podemos ilustrarlo para un sistema de $m=3$ ecuaciones en $n=2$ incógnitas.
$$\textbf{Combinación de columnas equivalentes a }b\qquad\begin{bmatrix}1 & 0 \\ 5 & 4 \\ 2 & 4 \end{bmatrix}\begin{bmatrix} u \\ v \end{bmatrix}=\begin{bmatrix}b_1\\ b_2\\ b_3\end{bmatrix}$$
Con $m>n$ tenemos mayor cantidad de ecuaciones que de incógnitas (sistema compatible indeterminado). Por lo tanto, el conjunto solución del sistema $A\,x=b$ puede ser expresado como una combinación de las columnas de $A$. De este modo tenemos que
$$\textbf{Combinación de columnas }\qquad u\begin{bmatrix} 1 \\ 5 \\ 2 \end{bmatrix}+v\begin{bmatrix} 0 \\ 4 \\ 4 \end{bmatrix}=\begin{bmatrix} b_1\\ b_2 \\ b_3 \end{bmatrix}$$
La representación de este conjunto solución sería el siguiente:
![[Pasted image 20240814001407.png]]
### El espacio nulo de una matriz
Existe un caso particular para el sistema $A\,x=b$ "dual" al primero (espacio columna). Si suponemos que $b=0$, tenemos que siempre $x=0$ es solución del sistema, pero pueden existir infinitas soluciones para este sistema (siempre que sean mayor las incógnitas que el números de ecuaciones, $n>m$).
>[!note] El espacio nulo de una matriz
>Conjunto de vectores $x$ tal que $A\,x=0$. Este espacio es denotado como $N(A)$ y su dimensión es definida por la dimensión de las columnas de $A$.

Este conjunto cumple las propiedades de un subespacio puesto que es cerrado para la suma y que la ecuación homogénea tiene al menos la solución trivial ($x=0$). Si $A$ tiene más columnas que filas (i.e., más incógnitas que ecuaciones, $n>m$), es probable que haya infinitas soluciones.
#### Ejemplo con una Matriz $A$ y su Espacio Nulo
Considera la matriz:
$$
A = \begin{pmatrix}
1 & 0 \\
5 & 4 \\
2 & 4
\end{pmatrix}
$$
La ecuación $Ax = 0$ con $x = \begin{pmatrix} u \\ v \end{pmatrix}$ resulta en el sistema de ecuaciones:
$$
\begin{aligned}
u &= 0, \\
v &= 0.
\end{aligned}
$$
En este caso, el espacio nulo de $A$ es trivial; solo contiene el vector $(0, 0)$. Esto ocurre porque las columnas de $A$ son linealmente independientes, lo que implica que $A$ no tiene otras soluciones no triviales para $Ax = 0$.
#### Espacio Nulo más Grande: Una Matriz $B$
Ahora, considera una matriz $B$ que es una extensión de $A$:
$$
B = \begin{pmatrix}
1 & 0 & 1 \\
5 & 4 & 9 \\
2 & 4 & 6
\end{pmatrix}
$$
En este caso, la tercera columna de $B$ es una combinación lineal de las dos primeras (específicamente, es la suma de las dos primeras columnas). Esto amplía el espacio nulo, ya que ahora hay más soluciones para $Bx = 0$. El espacio nulo de $B$ incluye todos los vectores de la forma $(c, c, -c)$, donde $c$ es cualquier número real. Este conjunto forma una línea que pasa por el origen en el espacio tridimensional.
## Sistemas generadores e independencia lineal
### Combinación lineal
>[!note] Combinación lineal
>Sea $V$ un $K$-espacio vectorial, y sea $G=\{v_1,...,v_r\}\subset V$. Una *combinación lineal* de $F$ es un elemento $v\in V$ tal que $v=\sum_{i=1}^r\alpha_i\cdot v_i$ con $\alpha_i\in K$ para cada $1\leq i\leq r$.

Podemos extender este concepto a los subconjuntos no necesariamente finitos del espacio vectorial considerado.

>[!note] Combinación lineal en subconjuntos no finitos
>Sea $V$ un $K$-espacio vectorial, sea $I$ un conjunto de índices y sea $G=\{v_i/i\in I\}\subset V$. Una *combinación lineal* de $G$ es un elemento $v\in V$ tal que $v=\sum_{i\in I}\alpha_i\cdot v_i$ donde $\alpha_i=0$ salvo para finitos $i\in I$. 
### Sistemas generadores
>[!note] Sistemas generadores
>Sea $V$ un $K$-espacio vectorial y sea $G\subset V$. Se dice que $F$ es un sistema de generadores de $V$ ( y se nota $<G> =V$) si todo elemento de $V$ es una combinación lineal de $G$.

Para cualquier sistema de generadores $G$ de un $K$-espacio vectorial $V$ que cumpla cierta propiedad adicional, que llamaremos independencia lineal, la cantidad de elementos de $G$ estará fija.
==Al considerar un sistema de generadores de un $K$-espacio vectorial $V$, hay que hallar sistemas de generadores que sean minimales respecto de la inclusión==, o sea, tal que ningún subconjunto propio sea también un sistema de generadores de $V$.
A continuación, se presentan resultados que caracterizan a los conjuntos con esta propiedad.
#### Resultado 1: Conjunto de vectores contenido en un subespacio específico
>[!note] Conjunto de vectores contenido en un subespacio específico
>Sean $V$ un $K$-espacio vectorial, $S$ un subespacio de $V$ y $\{v_1,...,v_n\}\subset V$. Entonces $<v_1,...,v_n>\subset S\iff v_i\in S,\forall 1\leq i\leq n$.

Esto establece una relación entre un conjunto de vectores y un subespacio vectorial.
Dado un espacio vectorial $V$ sobre un campo $K$, un subespacio $S$ de $V$, y un conjunto de vectores $\{v_1,...,v_n\}$ que pertenecen a $V$, el subespacio generado por estos vectores está contenido en $S$ si y sólo si cada vector $v_i$ está en $S$ para todo $i$ en el rango de $1$ a $n$. Por lo tanto:
- Si cada vector individual $v_i$ está en $S$, entonces cualquier combinación lineal de estos vectores también estará en $S$. Esto es debido a que los subespacios vectoriales son cerrados bajo la suma y la multiplicación por escalares.
- Si el subespacio generado por $\{v_1,...,v_n\}$ está contenido en $S$, entonces cada vector individual $v_i$ debe estar en $S$. Esto se debe a que cada $v_i$ es una combinación lineal de los vectores en $\{v_1,...,v_n\}$.
>[!tip] Verificación de conjunto contenido en un subespacio específico.
>Esta proposición nos da una forma de verificar si un conjunto de vectores está contenido en un subespacio específico.
>- Si puedes confirmar que cada vector individual está en el subespacio, entonces sabes que todas las posibles combinaciones lineales de esos vectores también estarán en el subespacio.
>- Y viceversa, si sabes que todas las combinaciones lineales están en el subespacio, entonces puedes asegurar que cada vector individual también está en el subespacio.
##### Demostración
###### $\implies$
Para cada $1\leq i\leq n$,$$v_i=0\cdot v_1+...+0\cdot v_{i-1}+1\cdot v_i+0\cdot v_{i+1}+...+0\cdot v_n\in<v_1,...,v_n>\subset S,$$de donde $v_i\in S$.
###### $\Longleftarrow$
Como $v_1,...,v_n\in S$ y $S$ es un subespacio, entonces $\sum_{i=1}^n\alpha_i v_i\in S\forall\alpha_i\in K$. Luego, $<v_1,...,v_n>$.
#### Resultado 2: Relación entre generadores y subespacio generado
>[!note] Relación generadores-subespacio generado
>Sea $V$ un $K$-espacio vectorial, y sea $\{v_1,...,v_n,v_{n+1}\}\subset V$. Entonces $<v_1,...,v_n,v_{n+1}>=<v_1,...,v_n>\iff v_{n+1}\in<v_1,...,v_n>$.

Dado un espacio vectorial $V$ sobre un campo $K$, y un conjunto de vectores $\{v_1,...,v_n,v_{n+1}\}$ que pertenecen a $V$, el subespacio generado por estos es igual al espacio generado por $\{v_1,...,v_n\}$ si y sólo si el vector $v_{n+1}$ está en el subespacio generado por $\{v_1,...,v_n\}$.
Entonces,
- Si el vector $v_{n+1}$ puede ser expresado como una combinación lineal de los vectores $\{v_1,...,v_n\}$ (o sea, $v_{n+1}\in\{v_1,...,v_n\}$), entonces agregar $v_{n+1}$ al conjunto de generadores no cambia el subespacio generado.
- Si el subespacio generado por $\{v_1,...,v_n,v_{n+1}\}$ es igual al subespacio generado por $\{v_1,...,v_n\}$, entonces, el vector $v_{n+1}$ debe estar en el subespacio generado por $\{v_1,...,v_n\}$. Esto se debe a que si $v_{n+1}$ no estuviera en $<v_1,...,v_n>$, entonces agregarlo al conjunto de generadores ampliaría el subespacio generado.
>[!tip] Verificación de $v_{n+1}\in<v_1,...,v_n>$
>Esta proposición nos permite verificar si un vector está en el subespacio generado por un conjunto de vectores.
>- Si puedo confirmar que dicho vector está en el subespacio, entonces agregarlo al conjunto no cambiará el subespacio generado.
>- Si sé que agregar el vector al conjunto de generadores no cambia el subespacio generado, entonces puedo estar seguro de que el vector está en el subespacio.
##### Demostración
###### $\implies$
Se tiene $<v_1,...,v_n,v_{n+1}>\subset<v_1,...,v_n>$. Entonces, por la proposición anterior, $v_{n+1}\in<v_1,...,v_n,v_{n+1}>$.
###### $\Longleftarrow$
Por hipótesis, $v_{n+1}\in<v_1,...,v_n>$. Además $v_i\in<v_1,...,v_i>\forall1\leq i\leq n$. Entonces, $<v_1,...,v_n,v_{n+1}>\subset<v_1,...,v_n>$. Por otro lado, $v_i\in<v_1,...,v_n,v_{n+1}>\forall1\leq i\leq n$, entonces vale $$<v_1,...,v_n,v_{n+1}>\supset<v_1,...,v_n>$$Luego $<v_1,...,v_n,v_{n+1}>=<v_1,...,v_n>$.
### Independencia lineal
>[!note] Independencia lineal
>Sea $V$ un $K$-espacio vectorial y sea $\{v_\alpha\}_{\alpha\in I}$ una familia de vectores de $V$. Se dice que $\{v_\alpha\}_{\alpha\in I}$ es **linealmente independiente** (l.i.) si
>$$\sum_{\alpha\in I}a_\alpha\cdot v_\alpha=0\implies a_\alpha=0,\forall \alpha\in I$$
>Caso contrario, $\{v_\alpha\}_{\alpha\in I}$ no es linealmente independiente, se dice que es **linealmente dependiente** (l.d.).

%%
Aunque, a diferencia de un conjunto, una familia puede contener elementos repetidos, utilizaremos ambos términos indistintamente, aceptando repeticiones.
%%
#### Criterio de l.i. 1
La noción de independencia lineal está íntimamente relacionada con la minimalidad de un sistema de generadores. Más precisamente:
>[!note] Condiciones de un conjunto linealmente independiente
>Sea $V$ un $K$-espacio vectorial y sean $v_1,...,v_n\in V$. Entonces el conjunto $\{v_1,...,v_n\}$ es linealmente independiente si y sólo si
>$$<v_1,...,v_n>\neq<v_1,...,\hat{v_i},...,v_n>\ \forall1\leq i\leq n$$
##### Demostración
###### $\implies$
Supongamos que $<v_1,...,\hat{v_i},...,v_n>=<v_1,...,v_n>$. En particular: $v_i\in<v_1,...,\hat{v_i},...,v_n>$, es decir, existen $\alpha_j\in K(j\neq i)$ tales que $v_i=\sum_{j=1, j\neq i}^n\alpha_jv_j$. Entonces: $$0=\sum_{j=1}^{i-1}\alpha_jv_j+(-1)v_i+\sum_{j=i+1}^n\alpha_jv_j,$$de donde $\{v_1,...,v_n\}$ no es linealmente independiente.
###### $\Longleftarrow$
Si $\{v_1,...,v_n\}$ es linealmente dependiente, existen $\alpha_1,...,\alpha_n\in K$ no todos nulos, tales que $\sum_{i=1}^n\alpha_iv_i=0$. Sin pérdida de generalidad, supongamos que $\alpha_n\neq0$. Entonces, $$v_n=-\sum_{i=1}^{n-1}\frac{\alpha_i}{\alpha_n}\cdot v_i\in<v_1,...,v_{n-1}>$$Luego, $<v_1,...,v_n>=<v_1,...,v_n>$.
#### Criterio de l.i. 2
>[!note] Propiedades invariantes de la independencia lineal bajo operaciones elementales
>Sea $V$ un $K$-espacio vectorial. Entonces:
>1. $\{v_1,...,v_i,...,v_j,...,v_n\}\subset V$ es l.i. $\iff \{v_1,...,v_j,...,v_i,...,v_n\}\subset V$ es l.i
>2. $\{v_1,...,v_i,...,v_n\}\subset V$ es l.i. $\iff \{v_1,...,\lambda v_i,...,v_n\}\subset V$ es l.i. para $\lambda\in K-{0}$.
>3. $\{v_1,...,v_i,...,v_j,...,v_n\}\subset V$ es l.i. $\iff\{v_1,...,v_i+\lambda v_j,...,v_j,...,v_n\}\subset V$ es l.i. para $\lambda\in K$.
##### Demostración
1. Se deduce del hecho que en un conjunto no interesa el orden de sus elementos.
2. Supongamos que $\{v_1,...,v_i,...,v_n\}$ es l.i. Sean $\alpha_1,...,\alpha_n\in K$ tales que $\alpha_1v_1+...+\alpha(\lambda v_i)+...+\alpha_nv_n=0$. Entonces se tiene que $\alpha_j=0$ para cada $j\neq i$ y que $\alpha_i\cdot\lambda=0$. Puesto que $\lambda\neq0$, resulta que también $\alpha_i=0$. Luego, el conjunto $\{v_1,...,\lambda v_i,...,v_n\}$ es linealmente independiente. Esto prueba la equivalencia, puesto que para demostrar la otra implicación basta multiplicar el i-ésimo vector del conjunto por $\frac{1}{\lambda}$.
3. Supongamos que $\{v_1,...,v_i,...,v_j,...,v_n\}$ es linealmente independiente. Sean $\alpha_1,...,\alpha_n\in K$ tales que:$$0=\alpha_1v_1+...+\alpha_i(v_i+\lambda v_j)+...+\alpha_jv_j+...+\alpha_nv_n$$$$=\alpha_1v_1+...+\alpha_iv_i+...+(\alpha_i\lambda+\alpha_j)v_j+...+\alpha_nv_n$$La independencia lineal de $\{v_1,...,v_i,...,v_j,...,v_n\}$ implica que $$\alpha_1=...=\alpha_i=...=\alpha_i\lambda+\alpha_j=...=\alpha_n=0$$de donde $\alpha_k=0$ para todo $1\leq k\leq n$. En consecuencia, el conjunto $\{v_1,...,v_i+\lambda v_j,...,v_j,...,v_n\}$ es l.i. La otra implicación se deduce de ésta, observando que el conjunto $\{v_1,...,v_n\}$ se obtiene de $\{v_1,...,v_i+\lambda v_j,...,v_j,...,v_n\}$ cambiando el $i$-ésimo vector $v_i+\lambda v_j$ por $(v_i+\lambda v_j)+(-\lambda)v_j=v_i$.
#### Criterio 3 (consecuencia del 2)
>[!note] Determinar l.i. mediante la matriz de las generatrices.
>Para decidir si un subconjunto de vectores $\{v_1,...,v_r\}$ de $K^n$ es l.i. podemos proceder:
>- Considerar la matriz $A$ cuyas filas son los vectores $v_1,...,v_r$.
>- Triangular la matriz $A$.
>**Si la matriz obtenida tiene alguna fila nula, el conjunto es l.d.**, de lo contrario, es l.i.
### Bases
>[!note] Base de un espacio vectorial
>Sea $V$ un $K$-espacio vectorial. Una familia $\{v_\alpha\}_{\alpha\in I}$ es **base** del espacio vectorial $V$ si $\{v_\alpha\}_{\alpha\in I}$  es una familia linealmente independiente de $V$ que satisface $$<v_\alpha>_{\alpha\in I}=V$$
#### Ejemplos
1. En $K^n,B=\{e_1,...,e_n\}$, donde $(e_i)_i=1$ y $(e_i)_j=0$ si $j\neq i$, es una base llamada **base canónica de** $K^n$.
2. En $K^{n\times m},B=\{E^{ij}/1\leq i\leq n, 1\leq j\leq m\}$ es una base.
3. En $K[X],B=\{X^i/i\in\mathbb{N}_0\}$ es una base
#### Lema de intercambio
>[!note] Lema de intercambio
>Si $\mathcal{B}=\{v_1,...,v_n\}$ es una base de $V$ y $w\in V$, entonces existe $1\leq i\leq n$ tal que reemplazando en $\mathcal{B}$ a $v_i$ por $w$, el conjunto resultante sigue siendo una base de $V$.
### Dimensión
#### Teorema de Dimensión de una Base en un espacio vectorial
>[!note] Teorema de Dimensión de una Base en un espacio vectorial
>Sea $V$ un $K$-espacio vectorial. Supongamos que $<v_1,...,v_r> =V$ y que $\{w_1,...,w_s\}\subset V$ es una familia linealmente independiente. Entonces $s\leq r$.
##### Demostración
Como $V=<v_1,...,v_r>$, para cada $1\leq i\leq s$, existen $\alpha_{ij}\in K(\leq j\leq r)$ tales que $w_i=\sum_{j=1}^r\alpha_{ij}v_j$. Consideremos el siguiente sistema de $r$ ecuaciones y $s$ incógnitas:
$$\sum_{h=1}^s\alpha_{hj}x_h=0\quad 1\leq j\leq r$$
Sea $(\beta_1,...,\beta_s)$ una solución del sistema. Entonces:
$$\sum_{h=1}^s\beta_hw_h=\sum_{h=1}^s\beta_h\left(\sum_{j=1}^r\alpha_{hj}v_j\right)=\sum_{h=1}^s\left(\sum_{h=1}^r\beta_h\alpha_{hj}v_j\right)$$
$$=\sum_{j=1}^r\left(\sum_{h=1}^s\beta_h\alpha_{hj}v_j\right)=\sum_{j=1}^r\left(\sum_{h=1}^s\beta_h\alpha_{hj}\right)v_j=0$$
Dado que $\{w_1,...,w_s\}$ es l.i., debe ser $(\beta_1,...,\beta_s)=0$.
En consecuencia, el sistema $rs$ tiene solución única, de donde se deduce que la cantidad de ecuaciones del sistema es mayor o igual que el número de variables, es decir $r\geq s$.
#### Igualdad de dimensiones de bases de espacios vectoriales
>[!info] Igualdad de Dimensiones de Bases de dos espacios vectoriales
>Sea $V$ un $K$-espacio vectorial, y sean $B_1$ y $B_2$ dos bases de $V$. Si $B_1=\{w_1,...,w_n\}$ y $B_2=\{v_1,...,v_m\}$, entonces $n=m$.
##### Demostración
Por el teorema anterior:
- Sean $B_1$ un sistema de generadores de $V$ y $B_2$ un conjunto l.i., entonces $n\geq m$.
- Sean $B_2$ un sistema de generadores de $V$ y $B_1$ un conjunto l.i., entonces $m\geq n$.
Luego, $n=m$.
#### Definición de dimensión de un espacio vectorial
>[!note] Dimensión de un espacio vectorial
>Sea $V$ un $K$-espacio vectorial y sea $B=\{v_1,...,v_n\}$ una base de $V$. Diremos entonces que $n$ es la **dimensión** de $V$. En este caso, diremos que $V$ es un $K$-espacio vectorial de *dimensión finita*, para distinguirlo de los e.v. que no admiten una base con finitos elementos.
>Por convención, la dimensión de {0} es 0.
##### Notación de dimensión de un espacio vectorial
Si $n$ es la dimensión del $K$-espacio vectorial $V$, escribimos $n=\text{dim}_K V$, o simplemente $\text{dim}V$ si el cuerpo $K$ queda claro por el contexto.
#### Propiedad de Representación única de vectores en un e.v. 
>[!note] Representación vectorial única en una Base de un espacio vectorial
>Sea $V$ un $K$-espacio vectorial de dimensión finita. Sea $\{v_1,...,v_n\}$ una base de $V$. Entonces para cada $x\in V$ existen únicos $\alpha_1,...,\alpha_n\in K$ tales que $x=\sum_{i=1}^n\alpha_iv_i$.
##### Demostración
La existencia se deduce de que, por ser una base de $V$, $\{v_1,...,v_n\}$ es un sistema de generadores de $V$.
Supongamos que $\sum_{i=1}^n\alpha_iv_i=\sum_{i=1}^n\beta_iv_i$, entonces $\sum_{i=1}^n(\alpha_i-\beta_i)v_i=0$. Como $\{v_1,...,v_n\}$ es un conjunto l.i., $\alpha_i-\beta_i=0,\forall 1\leq i\leq n$. Luego, $\alpha_i=\beta_i,\forall 1\leq i\leq n$, lo que prueba la unicidad.
### Cálculo de bases de e.v. de dimensiones finitas
>[!note] Cálculo de base de un espacio vectorial de dimensión finita
>Sea $V$ un $K$-espacio vectorial de dimensión finita.
>1. Sea $\{v_1,...,v_s\}$ un sistema de generadores de $V$. Entonces existe un subconjunto $G\subset\{v_1,...,v_s\}$ que es una base de $V$.
>2. Sea $\{w_1,...,w_r\}$ un conjunto linealmente independiente de $V$. Entonces existen elementos $w_{r+1},...,w_n\in V$ tales que $\{w_1,...,w_r,w_{r+1},...,w_n\}$ es una base de $V$.
#### Demostración
1. Si $\{v_1,...,v_s\}$ es l.i., entonces es una base de $V$. Si no es l.i., alguno de los vectores del conjunto es combinación lineal de los otros. Supongamos que $v_s\in<v_1,...,v_{s-1}>$. Consideramos ahora $\{v_1,...,v_{s-1}\}$, que es un sistema de generadores de $V$, y procedemos inductivamente.
2. Sea $B=\{z_1,...,z_n\}$ una base de $V$. Sea $G_0=<w_1,...,w_r>$. Consideramos $$G_1:=\begin{cases}\{w_1,...,w_r,z_1\}\ \text{si}\ z_1\not\in<G_0>\\\{w_1,...,w_r\}\ \text{si}\ z_1\in<G_0>\end{cases}$$Se produce inductivamente para $2\leq i\leq n$, es decir: $$G_i:=\begin{cases}G_{i-1}\cup\{z_i\} \ \text{si}\ z_1\not\in<G_0>\\ G_{i-1}\ \text{si}\ z_i\in<G_{i-1}>\end{cases}$$Observar que $\{w_1,...,w_r\}\subset G_i,\forall 1\leq i\leq n$. Además para cada $1\leq i\leq n$, $<z_1,...,z_i>\subset<G_i>$, y $G_i$ es un conjunto l.i. En particular, $V=<z_1,...,z_n>\subset<G_n>$ y $G_n$ es l.i. Luego, $G_n$ es una base de $V$.
## Suma y suma directa
### Suma de subespacios
Sean $U$ y $W$ subconjuntos de un e.v. $V$. La suma de $U$ y $W$, escrito $U+W$, consiste en todas las sumas $u+w$ en las que $u\in U$ y $w\in W$. Esto es,
$$U+W=\{u+w:u\in U, w\in W\}$$
#### Suma cerrada
>[!note] Suma cerrada para subespacios
>Dada la suma $U+W$ de los subespacios $U$ y $W$ de $V$ es también un subespacio de $V$.

Recordemos que la intersección $U\cap W$ es subespacio de $V$.
#### Teorema de la Dimensión de la Suma de Subespacios
>[!note] Teorema de la Dimensión de la Suma de Subespacios
>Sean $U$ y $W$ subespacios de dimensión finita de un e.v. $V$. En ese caso, $U+W$ tiene dimensión finita y
>$$\dim(U+W)=\dim U+\dim W-\dim(U\cap W)$$
### Suma directa
>[!note] Suma directa
>Se dice que el e.v. $V$ es la **suma directa** de sus subespacios $U$ y $W$, denotado por $$V=U\oplus W$$si todo vector $v\in V$ puede escribirse de una y sólo una forma como $v=u+w$, con $u\in U$ y $w\in W$.

De este modo, un e.v. $V$ es la suma directa de sus subespacios $U$ y $W$ si y sólo si:
- $V=U+W$
- $U\cap W=\{0\}$
#### Sumas directas generales
Podemos extender la noción de una suma directa a más de un sumando. Esto es $V$ es la suma directa de los subespacios $W_1,W_2,...,W_r$, escrito
$$V=W_1\oplus W_2\oplus\cdots\oplus W_r$$
Si todo vector $v\in V$ puede escribirse de una y sólo una forma como
$$v=w_1+w_2+\cdots+w_r$$
donde $w_1\in W_1,w_2\in W_2,...,w_r\in W_r$.
#### Teorema de la Independencia y Dimensión en Sumas Directas
>[!note] Teorema de la Independencia y Dimensión en Sumas Directas
>Sea $V=W_1\oplus W_2\oplus\cdots\oplus W_r$. Además, para cada $i$, supongamos que $S_i$ es un subconjunto linealmente independiente de $W_i$. En tal caso,
>1. La unión $S=\cup_iS_i$ es linealmente independiente en $V$.
>2. Si $S_i$ es una base de $W_i$, entonces $S=\cup_iS_i$ es una base de $V$.
>3. $\dim V=\dim W_1+\dim W_2+\cdots+\dim W_r$.
#### Teorema de la Igualdad de Dimensiones y Suma Directa
>[!note] Teorema de la Igualdad de Dimensiones y Suma Directa
>Supongamos que $V=W_1+W_2+\cdots+W_r$ (donde $V$ tiene dimensión finita) y que $$\dim V=\dim W_1+\dim W_2+\cdots+\dim W_r$$Entonces $V=W_1\oplus W_2\oplus\cdots\oplus W_r$.
# Órdenes de Convergencia
## Secuencias convergentes
[[Clase 5 - Límite. Continuidad]]
>[!note] Secuencias convergentes
>Una secuencia $\{x_n\}$ converge a un número real $L$ si:
>$$\lim_{n\to\infty}x_n=L$$
>es decir, para todo $\epsilon>0$ existe un número real $r$ tal que:
>$$|x_n-L|<\epsilon,\forall n>r$$

>[!tip] Importancia de las secuencias convergentes
>Es importante entender la convergencia de secuencias numéricas puesto que hay aplicaciones prácticas que incluyen métodos iterativos para resolver sistemas de ecuaciones lineales o hallar valores propios de matrices, donde la convergencia garantiza la precisión numérica de los resultados computacionales. 
## Complejidad algorítmica
[[Teórica 08 - Complejidad de Algoritmos]]
## Órdenes de Convergencia
En álgebra lineal y análisis numérico, se utiliza una terminología especial para describir la rapidez con la que una secuencia converge. Sea $\{x_n\}$ una secuencia de números reales que tiende a un límite $x^*$. Se dice que la tasa de convergencia es al menos **lineal** si existen una constante $c < 1$ y un entero $N$ tales que:
$$
|x_{n+1} - x^*| \leq c |x_n - x^*| \quad \text{(para } n \geq N\text{)}
$$
Se dice que la tasa de convergencia es al menos **superlineal** si existe una secuencia $\epsilon_n$ tendiendo a 0 y un entero $N$ tales que:
$$
|x_{n+1} - x^*| \leq \epsilon_n |x_n - x^*| \quad \text{(para } n \geq N\text{)}
$$
La convergencia es al menos **cuadrática** si existen una constante $C$ (no necesariamente menor que 1) y un entero $N$ tales que:
$$
|x_{n+1} - x^*| \leq C |x_n - x^*|^2 \quad \text{(para } n \geq N\text{)}
$$
En general, si existen constantes $C$ y $\alpha$ y un entero $N$ tales que:
$$
|x_{n+1} - x^*| \leq C |x_n - x^*|^\alpha \quad \text{(para } n \geq N\text{)}
$$
se dice que la tasa de convergencia es de **orden $\alpha$**, al menos.
## Teorema de la función implícita
![[Clase 10 - Dirección de máximo crecimiento. Curvas de nivel revisadas. Teorema de la Función Implícita#Teorema de la Función Implícita|Teorema de la Función Implícita]]
![[Clase 11 - Teorema de la Función implícita en R3, gradiente como vector normal a superficie de nivel, etc.#Teorema de la Función implícita en $ mathbb{R} 3$|Teorema de la Función implícita en R3]]