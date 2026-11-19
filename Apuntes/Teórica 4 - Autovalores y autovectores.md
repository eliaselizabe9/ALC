# Autovalores y autovectores
## Definición de autovalor y autovector
>[!note] Autovector y autovalor
>Dada una matriz cuadrada $A\in\mathbb{K}^{n\times n}$, un vector $v\in\mathbb{K}^n$ no nulo tal que
>$$Av=\lambda v,\ \text{para algún}\ \lambda\in\mathbb{K}$$
>se llama **autovector** de $A$ con **autovalor** $\lambda$.
### Ejemplo
Buscamos los autovalores y autovectores de la matriz $A=\begin{pmatrix}2 & 3 \\ 2 & 1\end{pmatrix}$. Tenemos que encontrar $v$ y $\lambda$ tales que $Av=\lambda v$. El término de la derecha podemos reescribirlo como $\lambda Iv$, y esto nos permite agrupar:
$$(A-\lambda I)v=0$$
que en nuestro caso queda
$$\begin{pmatrix}2-\lambda & 3 \\ 2 & 1-\lambda \end{pmatrix}\begin{pmatrix}v_1\\v_2\end{pmatrix}=\begin{pmatrix}0\\0\end{pmatrix}$$
Si el determinante de $A-\lambda I$ es no nulo, el sistema tiene solución única $(v_1,v_2)=(0,0)$. Como estamos buscando vectores no nulos, debemos encontrar $\lambda$ tal que $\text{det}(A-\lambda I)=0$.
Calculamos
$$\text{det}(A-\lambda I)=(2-\lambda)(1-\lambda)-6=\lambda^2-3\lambda-4$$
Esta ecuación es un polinomio de grado 2 en $\lambda$ y tiene raíces
$$\lambda_1=-1\qquad \text{y}\qquad \lambda_2=4$$
Estos son los autovalores de $A$.
Para calcular los autovectores correspondientes a esos autovalores, calculamos los núcleos de las matrices
$$A-(-1)I_2=\begin{pmatrix}3 & 3\\ 2 & 2\end{pmatrix}\quad\text{y}\quad A-(4)I_2=\begin{pmatrix}-2 & 3 \\ 2 & -3\end{pmatrix}$$

``` Python
import numpy as np
import scipy . linalg

A = np.array([[2,3],[2,1]])
Al = A - (-1) * np.eye(2)
print(scipy.linalg.null_space(Al))
```
%%
\[\[-0,70710678]
\[0,70710678]]
%%
El comando `null_space` normaliza los vectores a norma-2 igual a 1. Podemos verificar que el espacio de autovectores correspondiente a $\lambda_1=-1$ es
$$E_{\lambda_1}=<(1,-1)>$$
y el espacio de autovectores correspondiente a $\lambda_2=4$ es
$$E_{\lambda_2}=<(3,2)>$$
Es decir, por ejemplo, $(1,-1)$ es autovector de $A$ de autovalor $-1$, como también lo es cualquier múltiplo no nulo de ese vector. En este caso, obtuvimos que el espacio de autovectores de cada autovalor tiene dimensión 1, o abusando un poco la notación, que a cada autovalor le corresponde un único autovector.
Repasando,
- Los *autovalores* de una matriz $A\in\mathbb{K}^{n\times n}$ son los valores $\lambda\in\mathbb{C}$ para los cuales $\text{det}(A-\lambda I_n)=0$.
- El espacio de autovalores de un autovector $\lambda$ es el núcleo de la matriz $A-\lambda I_n$
- Si $\lambda$ es autovalor, el núcleo siempre tiene dimensión mayor o igual a 1
- Si el núcleo tiene dimensión exactamente 1, decimos informalmente que el autovalor $\lambda$ tiene un único autovector.
## Polinomio característico
>[!note] Polinomio característico
>Dada $A\in\mathbb{K}^{n\times n}$, llamamos **polinomio característico** de $A$ al determinante
>$$x_A(\lambda)=\text{det}(\lambda I-A)$$

Los autovalores de $A$ son las raíces de $x_a(\lambda)$. El polinomio $\text{det}(\lambda I-A)$ es un polinomio mónico con las mismas raíces que $\text{det}(\lambda I-A)$, y llamamos indistintamente **polinomio característico** a cualquier de estos dos polinomios. 
### Ejemplo
Verificamos en `Python` los cálculos para las dos matrices $A$ y $B$ utilizando el comando `linalg.eig` de `numpy`.
``` Python
A = np.array([[2,3],[2,1]])
e = np.linalg.eig(A)   # e es un a lista con dos elementos
print("Autovalores: ", e[0])   # El primer elemento es un array de autovalores
print("Autovectores:\n", e[1]) # El segundo elemento es una matriz con los autovectores como columnas
```
%%
Autovalores: \[4. -1.]
Autovectores:
	\[\[0.83205029 -0.70710678]
	\[0.5547002 0.70710678]]
%%
### Invariancia de los autovalores bajo transposición
>[!note] Invariancia de los autovalores bajo transposición
>Para $A\in\mathbb{R}^{n\times n}$, los autovalores de $A$ y $A^T$ son los mismos.
#### Demostración
Recordemos que para cualquier $A\in\mathbb{R}^{n\times n}$, $\text{det}(A)=\text{det}(A^T)$
Los autovalores de $A$ son las raíces del polinomio característico $p_A(\lambda)=\text{det}(A-\lambda I)$ y los autovalores de $A^T$ son las raíces del polinomio característico $p_{A^T}(\lambda)=\text{det}(A^T-\lambda I)$. Como $(A-\lambda I)^T=A^T-\lambda I$, ambos determinantes son iguales y por lo tanto las raíces son las mismas.
### Autovalores nulos
>[!note] $\lambda =0$
>$\lambda=0$ es autovalor de $A\iff A$ no es inversible.
#### $\Longleftarrow$
Suponiendo que $A$ no es inversible. Entonces $\exists x\neq 0:Ax=0=0\cdot x$. Por lo tanto, $0$ es autovalor.
#### $\implies$
Si $\lambda=0$ es autovalor, entonces $\mathcal X_A(0)\implies\det(A-0I)=0\iff\det(A)=0$. Por lo tanto, $A$ no es inversible.
# Diagonalización
Podemos pensar una matriz $A\in\mathbb{K}^{n\times n}$ como una transformación lineal de $T_A:\mathbb{K}^n\to\mathbb{K}^n$ dada por $T_A(v)=A\cdot v$ para $v\in\mathbb{K}^n$. En este caso estamos pensando que tanto el vector $v$ como $T_A(v)$ están expresados en coordenadas de la base canónica. Veamos como cambia la matriz de una transformación cuando expresamos los vectores en otra base.
## Diagonalización de una matriz
>[!note] Diagonalización de una matriz
>Sea $A\in\mathbb{K}^{n\times n}$ tal que existe una base $\mathcal{B}=\{w_1,...,w_n\}$ de $\mathbb{K}^n$ formada por autovectores de $A$, con autovalores $\lambda_1,...,\lambda_n\in\mathbb{C}$ (no necesariamente distintos). Sea $f:\mathbb{K}^n\to\mathbb{K}^n$ la transformación asociada, $f(v)=Av$. La matriz de la $f$ en la base $\mathcal{B}$ es diagonal,
>$$[f]_B=\begin{pmatrix}\lambda_1 & & 0\\ & \ddots & \\ 0 & & \lambda_n\end{pmatrix}$$
### Ejemplo
Continuando el ejemplo anterior, consideramos $A=\begin{pmatrix}2 & 3 \\ 2 & 1\end{pmatrix}$ y la transformación asociada
$$f(v)=T_A(v)=Av=\begin{pmatrix}2 & 3 \\ 2 & 1\end{pmatrix}\begin{pmatrix}v_1\\v_2\end{pmatrix}$$
Para esta transformación vale
$$T((1,0))=(2,2)\qquad T((0,1))=(3,1)$$
También vimos que $w_1=(1,-1)$ y $w_2=(3,2)$ son autovectores de $A$ con autovalores $-1$ y $4$ respectivamente. Si tomamos la base $\mathcal{B}=\{w_1,w_2\}$, la transformación en esta base toma una forma más simple:
$$T(w_1)=-w_1\qquad T(w_2)=4w_2$$
o si escribimos a estos vectores en las coordenadas de la base $\mathcal{B}$, obtenemos
$$T((1,0)_\mathcal{B})=(-1,0)_\mathcal{B}\qquad T((0,1)_\mathcal{B})=(0,4)_\mathcal{B}$$
Podemos entonces considerar la matriz de la transformación $T$ expresada en coordenadas de la base $\mathcal{B}$ y llamamos $[f]_\mathcal{B}$ a esta matriz. Obtuvimos la diagonalización de la $f$ en la base $\mathcal{B}$
$$[f]_\mathcal{B}=\begin{pmatrix}-1 & 0 \\ 0 & 4\end{pmatrix}$$
### Relaciones con la matriz original: matrices semejantes
Sabiendo que dada una base $\mathcal{B}=\{w_1,...,w_n\}$ de $\mathbb{K}^n$ y un vector $v\in\mathbb{K}^n$ con coordenadas en la base $\mathcal{B}$
$$v=(u_1,...,u_n)_\mathcal{B}$$
podemos encontrar las coordenadas de $v$ en la base canónica mediante la multiplicación
$$v_\mathcal{E}=C_{\mathcal{B}\mathcal{E}}v_\mathcal{B}$$
donde $C_{\mathcal{B}\mathcal{E}}$ es la matriz de cambio de base de $\mathcal{B}$ a $\mathcal{E}$,
$$C_{\mathcal{B}\mathcal{E}}=\begin{pmatrix}| & | & & |\\b_1 & b_2 & \cdots & b_n\\ | & | & & |\end{pmatrix}$$
Si en cambio, dadas las coordenadas canónicas de un vector $v=(v_1,...,v_n)_\mathcal{E}$, queremos encontrar las coordenadas en la base $\mathcal{B}$, realizamos la multiplicación
$$v_\mathcal{B}=C_{\mathcal{E}\mathcal{B}}v_\mathcal{E}$$
donde $C_{\mathcal{E}\mathcal{B}}$ es la matriz de cambio de base de $\mathcal{E}$ a $\mathcal{B}$, que podemos calcular invirtiendo la matriz $C_{\mathcal{B}\mathcal{E}}$
$$C_{\mathcal{E}\mathcal{B}}=C_{\mathcal{B}\mathcal{E}}^{-1}$$
Por lo tanto, si $f$ es un a t.l. dada en coordenadas de la base canónica por una matriz $A$, podemos obtener la matriz de la transformación en una base $\mathcal{B}$ (es decir que tanto $v$ como $T_A(v)$ están expresados en coordenadas de la base $\mathcal{B}$) aplicando el siguiente cambio de base a la matriz $A$:
$$A_\mathcal{B}=C_{\mathcal{E}\mathcal{B}}AC_{\mathcal{B}\mathcal{E}}$$
Veamos que tanto $A$ como $A_\mathcal{B}$ representan la misma transformación lineal, pero expresada en distintas bases. En este caso, decimos que las matrices son **semejantes**.
>[!note] Matrices semejantes
>Dos matrices $A,B$ son **semejantes** si existe una matriz inversible $C$ tal que
>$$B=C^{-1}AC$$
### Ejemplo
Continuando el ejemplo $A=\begin{pmatrix}2 & 3 \\ 2 & 1\end{pmatrix}$. Vimos que $A$ tiene autovalores $\lambda_1=-1$ y $\lambda_2=4$ con autovectores $v_1=(1,-1)$ y $v_2=(3,2)$ respectivamente. Estos autovectores forman una base de $\mathbb{R}^2$. Si expresamos la matriz $A$ en la base $B=\{v_1,v_2\}$ obtenemos
$$A_\mathcal{B}=\begin{pmatrix}-1 & 0 \\ 0 & 4\end{pmatrix}$$
(en las columnas de $A_\mathcal{B}$ ponemos las imágenes de los vectores de la base $\mathcal{B}$ expresados también en coordenadas de la base $\mathcal{B}$).
Para verificar esto en `Python`, utilizamos el siguiente código:
``` Python
A = np.array([[2,3],[2,1]])
e = np.linalg.eig(A)
C_BE=e[1]
C_EB=np.linalg.inv(C_BE)
print("A_B = \n", C_EB @ A @ C_BE)
```
%%
A_B = 
	\[\[ 4.00000000e+00 0.00000000e+00]
	\[-5.55111512e-17 -1.00000000e+00]]
%%
### Teorema de Diagonalización de matrices
>[!note] Teorema de Diagonalización de matrices
>Si $A\in\mathbb{K}^{n\times n}$ y $\mathcal{B}=\{b_1,...,b_n\}$ es una base de $\mathbb{K}^n$ formada por autovectores de $A$ con autovalores $\{\lambda_1,...,\lambda_n\}$, entonces tomando
>- $D$ la matriz diagonal de autovalores $D=\begin{pmatrix}\lambda_1 & & \\ & \ddots & \\ & & \lambda_n\end{pmatrix}$
>- $C$ la matriz de autovectores de $A,C=\begin{pmatrix}| & & |\\b_1 & \cdots & b_n\\ | & & |\end{pmatrix}$ obtenemos la diagonalización de la matriz $A$: $$A=CDC^{-1}\qquad\text{y}\qquad D=C^{-1}AC$$
## Multiplicidad algebraica y geométrica
>[!note] Multiplicidad aritmética y geométrica
>Dada una matriz $A\in\mathbb{K}^{n\times n}$ y un autovalor $\lambda$ de $A$, definimos:
>- La **multiplicidad aritmética** de $\lambda$ es la multiplicidad de $\lambda$ como raíz del polinomio característico $\mathcal{X}_A$
>- La **multiplicidad geométrica** de $\lambda$ es la dimensión del espacio de autovectores asociado a $\lambda$, o equivalentemente, la dimensión de $\text{Nu}(A-\lambda I)$.
### Desigualdad de multiplicadores de autovalores
>[!note] Desigualdad de multiplicadores de autovalores
>Dada una matriz $A\in\mathbb{K}^{n\times n}$ y un autovalor $\lambda$ de $A$, la multiplicidad geométrica es siempre menor o igual que la multiplicidad aritmética.
#### Demostración
Sea $s$ la multiplicidad geométrica de $\lambda$ y $\{v_1,...,v_s\}$ una base de autovectores de $E_\lambda$, el espacio de autovectores asociados a $\lambda$. Completamos la base de $E_\lambda$ a una base de $\mathbb{K}^n$:
$$\mathcal{B}=\{v_1,...,v_s,w_1,...,w_{n-s}\}$$
Consideramos la matriz de la transformación $f$ definida por $A$ en la base $\mathcal{B}$:
$$[f]_\mathcal{B}=\left(\begin{array}{ccc|ccc} \lambda & & & & & \\ & \ddots & & & B & \\ & & \lambda & & &\\ \hline  &  &  &  \\  & 0 &  & & D & \\ & & & & & & \end{array}\right)$$
donde el primer bloque es una matriz de $s\times s$.
Mediante la fórmula para el cálculo de determinantes desarrollando una columna, es fácil ver que
$$\mathcal{X}_{[f]_\mathcal{B}}(x)=\text{det}(xI_n-[f]_\mathcal{B})=(x-\lambda)^s\text{det}(xI_n-[f]_\mathcal{B})$$
Ya vimos que $[f]_\mathcal{B}$ y $[f]_\mathcal{E}$ tiene el mismo polinomio característico. Por lo tanto, la multiplicidad algebraica de $\lambda$ como autovalor de $A$ es al menos $s$, como queríamos ver.
### Independencia lineal de autovectores correspondientes a autovalores distintos
Para una matriz $A\in\mathbb{K}^{n\times n}$, el polinomio característico tiene grado $n$ y por lo tanto $n$ raíces contando sus multiplicidades. Es decir, que la suma de las multiplicidades aritméticas de todos los autovalores es siempre exactamente $n$. Más aún, autovectores correspondientes a autovalores distintos son linealmente independientes.
>[!note] Independencia lineal de autovectores correspondientes a autovalores distintos
>Dada $A\in\mathbb{K}^{n\times n}$, sean $v_1,...,v_s$ autovectores de $A$ correspondientes a autovalores distintos $\lambda_1,...,\lambda_s$. El conjunto $\{v_1,...,v_s\}$ es linealmente independiente.
#### Demostración
Consideramos una relación de dependencia lineal:
$$a_1v_1+a_2v_2+\cdots+a_sv_s=0$$
y supongamos que es la relación que involucra la menor cantidad posible de vectores, en particular todos los $a_i$ son no nulos.
Si multiplicamos la igualdad por $(A-\lambda_1 I_n)$ obtenemos
$$a_1(A-\lambda_1 I_n)v_1+a_2(A-\lambda_1 I_n)v_2+\cdots+a_s(A-\lambda_1 I_n)v_2=0$$
y simplificando,
$$0+a_2(\lambda_2-\lambda_1)v_2+\cdots+a_s(\lambda_s-\lambda_1)v_s=0$$
que es la relación con menos términos no nulos que la anterior, lo cual es un absurdo.
### Diagonalizabilidad
Como consecuencia de la **Desigualdad de multiplicadores de autovalores** y la **Independencia lineal de autovectores correspondientes a autovalores distintos**, tenemos la definición de diagonalizabilidad de una matriz.
>[!note] Diagonalizabilidad
> Una matriz $A$ es diagonalizable si para todo autovalor $\lambda$ de $A$, la multiplicidad aritmética y la multiplicidad geométrica coinciden.
## Propiedades de las matrices diagonalizables
Si $A=C^{-1}DC$,
1. $A^2=C^{-1}DCC^{-1}DC=C^{-1}DDC=C^{-1}D^2C$
2. Inductivamente, $A^k=C^{-1}DC...C^{-1}DC=C^{-1}D^kC$
3. Si $A$ es inversibles, entonces $A^{-1}=(C^{-1}DC)^{-1}=C^{-1}D^{-1}C$
4. Si $p(x)\in\mathbb{R}[x],p(A)=C^{-1}p(D)C$
En todos estos casos, como vimos anteriormente, la función aplicada a una matriz diagonal se puede calcular fácilmente aplicando la función a cada elemento de la diagonal.
## Aplicación
%% Página 99 pdf %%
# Procesos de Markov
## Ejemplo introductorio
En una ciudad, todos los adultos tienen teléfono celular de una de las tres compañías A, B o C. Cada mes, algunos clientes se cambian de compañía, según la siguiente fórmula:
$$\begin{pmatrix}a_{k+1}\\b_{k+1}\\c_{k+1}\end{pmatrix}=\begin{pmatrix}0.9 & 0.15 & 0.25\\ 0.075 & 0.8 & 0.25 \\ 0.025 & 0.05 & 0.5\end{pmatrix}\begin{pmatrix}a_k \\ b_k \\ c_k\end{pmatrix}$$
Las casillas de la matriz indican que porcentaje de los clientes de cada compañía se queda en la misma compañía o se cambian a otra compañía. Por ejemplo, $a_{k+1}=0.9a_k+0.15b_k+0.25c_k$, lo que significa que el $90\%$ de los clientes de la compañía $A$ en el mes $k$ siguen en la compañía $A$ en el mes $k+1$, el $15\%$ de los clientes de la compañía $B$ se pasan a $A$ y el $25\%$ de los clientes de la compañía $C$ se pasan a $A$.
Queremos ver cuántos clientes s hay en cada compañía con el paso del tiempo, y si esas cantidades converge a algún valor límite. Para esto, contestamos las siguientes preguntas utilizando `Python`.
Si inicialmente hay 1000 clientes en cada compañía, es decir $v_0=(a_0,b_0,c_0)=(1000,1000,1000)$,
1. Hallar $v_1(a_1,b_1,c_1)$ y $v_2=(a_2,b_2,c_2)$
2. Hallar $v_{10}$ y $v_{100}$
3. Hallar (si existe) el límite de la cantidad de clientes en cada compañía, es decir, $\lim_{k\to\infty}(a_k,b_k,c_k)$.
Primero determinamos $v_1=Av_0$ y $v_2=Av_1$.
``` Python
import numpy as np

A = np.array([[0.9,0.15,0.25],[0.075,0.8,0.25],[0.025,0.05,0.5]])

v0 = np.array([1000,1000,1000])
v1 = A @ v0
print("v1 = ", v1)
v2 = A @ v1
print("v2 = ", v2)
```
%%
v1 = \[1300. 1125. 575.]
v2 = \[1482.5 1141.25 376.25]
%%
En general, para calcular $v^{(k)}$, tenemos
$$v^{(k)}=Av^{(k-1)}=A(Av^{(k-2)})=\cdots=A(A(\cdots(Av^{(0)})\cdots))$$
y podemos calcular estos vectores mediante una iteración en `Python`:
``` Python
v = np.array([1000,1000,1000])
k = 10
for i in range(k):
	v = A @ v
print("v(",k, ") = ", v)

k = 100
for i in range(k):
	v = A @ v
print("v (", k, ") = ", v)
```
%%
v(10) = \[1844.0598284 965.482631 190.4575406]
v(100) = \[1875. 937.5 187.5]
%%
``` Python
v = np.array([1000,1000,1000])
k = 1000
for i in range(k):
	v = A @ v
print("v(",k, ") = ", v)

k = 10000
for i in range(k):
	v = A @ v
print("v (", k, ") = ", v)
```
%%
v(1000) = \[1875. 937.5 187.5]
v(10000) = \[1875. 937.5 187.5]
%%
Observamos que $\lim_{k\to\infty}v^{(k)}=(1875,937.5,187.5)$.
Para poder analizar el comportamiento de $v^{(k)}$ o verificar que para cualquier vector inicial, existe $v^{(\infty)}=\lim_{k\to\infty}v^{(k)}$, necesitamos estudiar los procesos de Markov.
## Matrices de Markov (estocásticas)
>[!note] Matrices de Markov
>Decimos que una matriz $A\in\mathbb{R}^{n\times n}$ es **de Markov** (o estocástica) si
>- Todas las casillas de $A$ son números reales no negativos
>- La suma de las casillas en cualquier columna de $A$ es siempre 1

En particular, estas dos condiciones implican que todas las coordenadas de la matriz deben ser números entre 0 y 1.
### Espectro de las matrices de Markov
Observamos que si $v^\infty$ es el límite de un proceso de Markov, debe cumplir $Av^\infty=v^\infty$, es decir $v^\infty$ es un autovector de $A$ de autovalor 1.
Estudiamos ahora entonces los autovalores y autovectores de las matrices estocásticas.
#### Autovalor unitario de Matrices de Markov
>[!note] Autovalor Unitario de Matrices de Markov
>Si $A$ es una matriz de Markov, $\lambda=1$ es un autovalor de $A$.
##### Demostración
Como las columnas de $A$ suman 1, si $v=(1,1,1,...,1)$,
$$A^Tv=v$$
luego $v$ es autovector de $A^T$ con autovalor 1. Como los autovalores de $A$ y $A^T$ son los mismos concluimos que $\lambda=1$ es autovalor de $A$ (aunque no necesariamente el autovector correspondiente es $v=(1,1,...,1)$).
#### Cota del módulo de los autovalores de matrices de Markov
>[!note] Cota del módulo de los autovalores de matrices de Markov
>Si $A$ es una matriz de Markov y $\lambda$ es un autovalor de $A$, entonces $|\lambda|\leq1$.
##### Demostración
Si $Av=\lambda v$, para $1\leq i\leq n$
$$|\lambda||v_i|=|\lambda v_i|=|A_{i1}v_1+\cdots+A_{in}v_n|\leq|A_{i1}v_1|+\cdots+|A_{in}v_n|=A_{i1}|v_1|+\cdots+A_{in}|v_n|$$
(para la última igualdad utilizamos que todas las coordenadas de $A$ son no-negativas).
Ahora sumamos sobre todos los $i$:
$$\begin{array}{lll}
	|\lambda|(|v_1|+\cdots+|v_n|) & \leq & (A_{11}|v_1|+\cdots+A_{1n}|v_n|)+\cdots+(A_{n1}|v_1|+\cdots+A_{nn}|v_n|)\\
	 & = & (A_{11}+\cdots+A_{n1})|v_1|+\cdots+(A_{1n}+\cdots+A_{nn})|v_n|\\
	 & = & |v_1|+\cdots+|v_n|
\end{array}$$
Como $|v_1|+\cdots+|v_n|\geq0$, obtenemos $|\lambda|<1$.
#### Condición de suma cero para autovectores de autovalores no unitarios en matrices de Markov
>[!note] Condición de suma cero para autovectores de autovalores no unitarios en matrices de Markov
>Si $A$ es una matriz de Markov, $\lambda\neq1$ es un autovalor de $A$ y $v$ es un autovector asociado a $\lambda$, entonces $v_1+\cdots+v_n=0$.
##### Demostración
Si $Av=\lambda v$, para $1\leq i\leq n$,
$$\lambda v_i=A_{i1}v_1+\cdots+A_{in}v_n$$
y
$$\begin{array}{lll}
	\lambda(v_1+\cdots+v_n) & = & (A_{11}v_1+\cdots+A_{1n}v_n)+\cdots+(A_{n1}v_1+\cdots+A_{nn}v_n)\\
	 & = & (A_{11}+\cdots+A_{n1})v_1+\cdots+(A_{1n}+\cdots+A_{nn})v_n\\
	 & = & v_1+\cdots+v_n
\end{array}$$
Como $\lambda\neq1$, debe ser $v_1+\cdots+v_n=0$.
### Probabilidades y proporciones
Si consideramos un sistema con $n$ estados distintos, podemos considerar a las casillas $(i,j)$ de una matriz de Markov como la probabilidad de pasar del estado $j$ al estado $i$.
#### Ejemplo
Para modelar el estado del tiempo en una ciudad, consideramos tres estados: soleado, nublado y lluvioso, y la matriz de Markov que indica la probabilidad de pasar de un estado a otro de un día a otro:
$$\begin{pmatrix}
	0.6 & 0.2 & 0.1\\
	0.3 & 0.6 & 0.4\\
	0.1 & 0.2 & 0.5
\end{pmatrix}$$
Si consideramos que inicialmente tenemos la misma probabilidad de que el día sea soleado, nublado o lluvioso, es decir, comenzamos con un vector inicial de probabilidades:
$$v^{(0)}=\left(\frac{1}{3},\frac{1}{3},\frac{1}{3}\right)$$
calculando $v^{(k)}podremos estimar la probabilidad de que llueva en un día específico, o calculando $v^\infty$ podemos estimar la probabilidad de que llueva a largo plazo.
### Conjunto de estados
>[!note] Conjunto de estados
>Conjunto $\mathcal{P}_n$ de estados $$\mathcal{P}_n=\{\pi\in\mathbb{R}^n_+:\pi_1+\cdots+\pi_n=1\}$$que podemos pensar como vectores de probabilidades.
#### Invariancia del conjunto de estados bajo matrices de Markov
>[!note] Invariancia del conjunto de estados bajo matrices de Markov
>Si $\pi\in\mathcal{P}_n$ y $M\in\mathbb{R}^{n\times n}$ es una matriz de Markov, $M\pi\in\mathcal{P}_n$.

O sea, el conjunto de estados $\mathcal{P}_n$ es cerrado bajo la multiplicación por matrices de Markov, reflejando la estabilidad de la estructura probabilística.
#### Existencia de un vector estacionario para matrices de Markov
>[!note] Existencia de un vector estacionario para matrices de Markov
>Si $M\in\mathbb{R}^{n\times n}$, existe $\pi\in\mathcal{P}_n$ tal que $M\pi=\pi$.
##### Demostración (ejercicio)
Ver que si $\pi=(\pi_1,...,\pi_n)$ satisface $M\pi=\pi$, entonces también $|\pi|=(|\pi_1|,...,|\pi_n|)$ satisface $M|\pi|=|\pi|$.
### Estados de equilibrio
>[!note] Estado de equilibrio
>Un estado $\pi\in\mathcal{P}_n$, se dice **estado de equilibrio** de una matriz de Markov $M$ si cumple $$M\pi=\pi$$o equivalentemente, si es un autovector de $M$ de autovalor 1.
>- Toda matriz de Markov tiene al menos un estado de equilibrio.
>- Una matriz de Markov puede tener varios estados de equilibrio.
>- Para el autovalor $\lambda=1$, la la multiplicidad algebraica coincide con la multiplicad geométrica, y podemos hallar una base de autovectores de $E_1$ formada por vectores de estado.
### Estados límite
Nos preguntamos si dado un estado inicial $\pi^{(0)}$, existe
$$\pi^*=\lim_{k\to\infty}\pi^{(k)}$$
En tal caso,
$$\pi^*=\lim_{k\to\infty}\pi^{(k+1)}=\lim_{k\to\infty}M\pi^{(k)}=M\pi^*$$
es decir, $\pi^*$ es un estado de equilibrio.
Entonces, tenemos las siguientes preguntas:
- ¿Dado un estado inicial $\pi(0)$, el sistema evoluciona siempre hacia un estado límite?
- Si el estado de equilibrio no es único, ¿cómo determinamos a qué estado converge?
### Autovalores y autovectores
Supongamos que $M$ es diagonalizable, con autovalores $|\lambda_1|\geq|\lambda_2|\geq\cdots\geq|\lambda_n|\geq0$, y autovectores correspondientes $w_1,...,w_n$.
Usando las propiedades vistas, si $s$ es la multiplicidad del autovalor 1$, podemos suponer $\lambda_1=\lambda_2=\cdots=\lambda_s=1$ y $w_1,...,w_s$ vectores de estado.
Dado $\pi^{(0)}\in\mathcal{P}_n$, escribimos
$$\begin{array}{}
	\pi^{(0)}=a_1w_1+\cdots+a_nw_n\\
	\pi^{(k)}=a_1(\lambda_1)^kW_1+\cdots+a_n(\lambda_n)^kw_n
\end{array}$$
### Cálculo analítico de estados límite
A partir de la expresión
$$\begin{array}{}
	\pi^{(k)}=a_1(\lambda_1)^kW_1+\cdots+a_n(\lambda_n)^kw_n\\
	=a_1w_1+\cdots+a_sw_s+a_{s+1}(\lambda_{s+1})^kw_{s+1}+\cdots+a_n(\lambda_n)^kw_n
\end{array}$$
utilizando que $\lambda_i^k\to0$ si $|\lambda_i|<0$ deducimos
- Si s el único autovalor de módulo 1, $\lim_{k\to\infty}\pi^(k)=a_1w_1+\cdots+a_sw_s$
- Si además $s=1,\pi^{(k)}\to w_1$
- Si existe $\lambda_i\neq1,|\lambda_i|=1$, el proceso de Markov puede no tener límite (depende del vector inicial de estados).
### Diagonalización de matrices de Markov
Si $M$ es diagonalizable, con autovalores $1=|\lambda_1|\geq|\lambda_2|\geq\cdots\geq|\lambda_n|\geq0$, y autovectores correspondientes $w_1,...,w_n$, tenemos
$$\begin{array}{}
	M=C_{\mathcal{B}\mathcal{E}}DC_{\mathcal{E}\mathcal{B}}=C_{\mathcal{B}\mathcal{E}}\begin{pmatrix}
		\lambda_1 \\
		 & \ddots \\
		 & & \lambda_n
	\end{pmatrix}C_{\mathcal{E}\mathcal{B}}\\
	M^k=C_{\mathcal{B}\mathcal{E}}D^kC_{\mathcal{E}\mathcal{B}}=C_{\mathcal{B}\mathcal{E}}\begin{pmatrix}
		\lambda_1^k \\
		 & \ddots \\
		 & & \lambda_n^k
	\end{pmatrix}C_{\mathcal{E}\mathcal{B}}
\end{array}$$
con $C_{\mathcal{B}\mathcal{E}}=\begin{pmatrix}| & & | \\ w_1 & \cdots & w_n\\ | & & |\end{pmatrix}$.
### Matriz $M^\infty$
Notamos $M^\infty=\lim_{k\to\infty}M^k$, cuando ese límite existe.
#### Propiedades de las matrices $M^\infty$
Estas matrices cumplen las siguientes propiedades
- Si existe $M^\infty,\pi^\infty=M^\infty\pi^{(0)}$ para cualquier vector inicial de estados.
- Existe $M^\infty$ si y sólo si $\lambda=1$ es el único autovalor de módulo 1.
- En este caso, $M^\infty=C_{\mathcal{B}\mathcal{E}}D^\infty C_{\mathcal{E}\mathcal{B}}$, y tomando $1^\infty=1$ y $\lambda^\infty=0$ para $|\lambda|<1$ obtenemos
	$$D=\begin{pmatrix}
		\lambda_1^\infty \\
		 & \ddots \\
		 & & \lambda_n^\infty
	\end{pmatrix}=\begin{pmatrix}
		1 \\ & \ddots\\ & & 1\\ & & & 0
	\end{pmatrix}$$
- Si el autovalor $\lambda=1$ tiene multiplicidad 1, con autovector de estado $w$, entonces
	$$M^\infty=\begin{pmatrix}
		| & & |\\w & \cdots & w\\ | & & |
	\end{pmatrix}$$
# El método de la potencia (6)
>[!note] Método de la potencia
>Método para calcular mediante aproximaciones numéricas los autovalores de una matriz. Consiste en un algoritmo iterativo que permite encontrar el autovalor dominante y el autovector correspondiente de una matriz.
## Ejemplo introductorio
Consideremos la matriz
$$A=\begin{pmatrix}
	0.9 & 0.075 & 0.025 \\
	0.15 & 0.8 & 0.05 \\
	0.25 & 0.25 & 0.5
\end{pmatrix}$$
Tomamos un vector cualquiera no nulo, por ejemplo $v=(1,2,3)$ y calculamos $A^k v$ para distintos valores de $k$. O, equivalentemente, calculamos distintos términos de la sucesión definida recursivamente por
$$\begin{array}{}
	v^{(0)}=(1,2,3)\\
	v^{(k)}=Av^{(k-1)}
\end{array}$$
Lo programamos en `Pyhton`.
``` Python
def estado (A, v, k):
	for i in range(k):
		v = A @ v
	return (v)

A = np.array([[0.9, 0.075, 0.025], [0.15, 0.8, 0.05], [0.25, 0.25, 0.5]])
v = np.array([1,2,3])
print(estado(A, v, 1))
print(estado(A, v, 10))
print(estado(A, v, 100))
print(estado(A, v, 1000))
```
%%
\[1.125 1.9 2.25 ]
\[1.41762971 1.47373372 1.45503428]
\[1.4375 1.4375 1.4375]
\[1.4375 1.4375 1.4375]
%%
Observamos en el experimento que $v^{(k)}$ tiende a $w=(1.4375, 1.4375, 1.4375)$. Si esto es cierto, entonces $Aw=w$. Luego encontramos un autovector de autovalor 1. ¿Cómo podemos demostrar lo que observamos?
Calculamos autovalores y autovectores de $A$ y verificamos si $A$ es diagonalizable.
``` Python
e = np.linalg.eig(A)
print(e[0])
print(e[1])
```
%%
\[1. 0.74142136 0.45857864]
\[\[-0.57735027 -0.44371857 -0.03400257]
 \[-0.57735027 0.81130706 -0.13017638]
 \[-0.57735027 0.38065035 0.99090763]]
%%
Obtenemos que $A$ es diagonalizable con autovalores $\lambda_1=1$, $\lambda_2=0.74$, $\lambda_3=0.45$. Llamamos $\mathcal{B}=\{w_1,w_2,w_3\}$ a la base de autovectores.
Luego, podemos escribir a cualquier vector $v^{(0)}\in\mathbb{R}^3$ como combinación lineal de vectores de la base $\mathcal{B}$:
$$v^{(0)}=a_1w_1+a_2w_1+a_3+w_3$$
Usando esta representación, ¿quién es $Av^{(0)}$? Obtenemos
$$\begin{array}{lll}
	v^{(1)}=Av^{(0)} & = & a_1Aw_1+a_2Aw_2+a_3Aw_3\\
	 & = & a_1\lambda_1w_1+a_2\lambda_2w_2+a_3\lambda_3w_3
\end{array}$$
y si multiplicamos nuevamente por $A$
$$\begin{array}{lll}
	v^{(2)}=A^2v=A(Av) & = & a_1A\lambda_1w_1+a_2A\lambda_2w_2+a_3A\lambda_3w_3\\
	 & = & a_1\lambda_1^2w_1+a_2\lambda_2^2w_2+a_3\lambda_3^2w_3
\end{array}$$
Repitiendo el mismo razonamiento obtenemos
$$v^{(k)}=A^kv=a_1(\lambda_1)^kw_1+a_2(\lambda_2)^kw_2+a_3(\lambda_3)^kw_3$$
Como $|\lambda_2|<1$ y $|\lambda_3|<1$ tenemos que $(\lambda_2)^k\stackrel{k\to\infty}{\longrightarrow}0$ y $(\lambda_3)^k\stackrel{k\to\infty}{\longrightarrow}0$. Por otro lado, $\lambda_1^k=1,\forall k\in\mathbb{N}$.
Tomando límite, concluimos entonces
$$v^{(k)}\stackrel{k\to\infty}{\longrightarrow}a_1w_1$$
es decir $v^{(k)}$ tiende a un vector de autovalor 1, que es el autovalor de mayor módulo de la matriz.
En general, las matrices pueden tener autovalores de cualquier módulo y la sucesión anterior podría tender a 0 o divergir o no tener límite. Para el caso general, realizamos la siguiente modificación.
Tomamos un vector cualquiera no nulo $v^{(0)}$ y en cada iteración de la recursión anterior normalizamos el resultado a norma-2=1:
$$v^{(k)}=\frac{Av^{(k-1)}}{||Av^{(k-1)}||_2},\text{para }k\geq1$$
Calculamos en `Python` los vectores de esta sucesión tomando $v^{(0)}=(1,2,3)$ y 
$$A=\begin{pmatrix}
	3 & 1 & 2\\
	0 & 1 & -2\\
	1 & 2 & 4
\end{pmatrix}$$
``` Python
def estado (A, v, k):
	for i in range (k):
		Av = A @ v
		v = Av / np.linalg.norm
	return(v)

A = np.array([[3,1,2],[0,1,-2],[1,2,4]])
v = np.array([1,2,3])

print(estado(A, v, 1))
print(estado(A, v, 10))
print(estado(A, v, 100))
print(estado(A, v, 1000))
```
%%
\[ 0.53295174 -0.19380063 0.82365269]
\[ 0.74404261 -0.37049752 0.55599656]
\[ 0.74278135 -0.37139068 0.55708601]
\[ 0.74278135 -0.37139068 0.55708601]
%%
Si observamos que en este caso también existe límite $v^{(k)}\stackrel{k\to\infty}{\longrightarrow}w=(2.97,-1.49,2.23)$. ¿Será cierto en este caso también que $w$ es un autovector? ¿Cuál es el autovalor correspondiente? Para responder estas preguntas, utilizamos nuevamente que si $w$ es el límite de la sucesión, al aplicar la fórmula recursiva a $w$ obtenemos nuevamente $w$. Es decir,
$$w=\frac{Aw}{||Aw||_2}$$
y despejando, obtenemos
$$Aw=||Aw||_2w$$
Concluimos que $w$ es un autovector de $A$ con autovalor $\lambda=||Aw||_2=4$ (verificar con `Python`).
Para demostrar el resultado en este caso, calculamos nuevamente autovalores y autovectores de $A$:
``` Python
e = np.linalg.eig(A)
print(e[0])
print(e[1])
```
%%
\[4.+0.j 2.+1.j 2.-1.j]
\[\[ 0.74278135+0.j 0.35355339-0.35355339j 0.35355339+0.35355339j]
 \[-0.37139068+0.j -0.70710678+0.j -0.70710678-0.j ]
 \[ 0.55708601+0.j 0.35355339+0.35355339j 0.35355339-0.35355339j]]
%%
Nuevamente la matriz es diagonalizable con autovalores $\lambda_1=4,\lambda_2=2+i,\lambda_3=2-i$. Aplicando el razonamiento anterior,
$$\begin{array}{l}
	v^{(1)}=\frac{Av^{(0)}}{||Av^{(0)}||_2}\\
	v^{(2)}=\frac{Av^{(1)}}{||Av^{(1)}||_2}=\frac{A\frac{Av^{(0)}}{||Av^{(0)}||_2}}{\left|\left|A\frac{Av^{(0)}}{||Av^{(0)}||_2}\right|\right|_2}=\frac{\frac{A^2v^{(0)}}{||Av^{(0)}||_2}}{\frac{||A^2v^{(0)}||_2}{||Av^{(0)}||_2}}=\frac{A^2v^{(0)}}{||A^2v^{(0)}||_2}\\
	\cdots\\
	v^{(k)}=\frac{Av^{(k-1)}}{||Av^{(k-1)}||_2}=\frac{A^kv^{(0)}}{||A^kv^{(0)}||_2}
\end{array}$$
Obtenemos en este caso
$$\begin{array}{lll}
	v^{(k)} & = &\frac{A^kv^{(0)}}{||A^kv^{(0)}||_2}\\
	 & = & \frac{a_1(\lambda_1)^kw_1+a_2(\lambda_2)^kw_2+a_3(\lambda_3)^kw_3}{||a_1(\lambda_1)^kw_1+a_2(\lambda_2)^kw_2+a_3(\lambda_3)^kw_3||_2}\\
	 & = & \frac{\lambda_1^k\left(a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+a_3\left(\frac{\lambda_3}{\lambda_1}\right)^kw_3\right)}{|\lambda_1|^k\left|\left|a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+a_3\left(\frac{\lambda_3}{\lambda_1}\right)^kw_3\right|\right|_2}\\
	 & = & \left(\frac{\lambda_1}{|\lambda_1|}\right)^k\frac{a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+a_3\left(\frac{\lambda_3}{\lambda_1}\right)^kw_3}{\left|\left|a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+a_3\left(\frac{\lambda_3}{\lambda_1}\right)^kw_3\right|\right|_2}
\end{array}$$
En nuestro ejemplo, $\lambda_1=4$ y $\left(\frac{\lambda_1}{|\lambda_1|}\right)^k=1$. Los cocientes $\frac{\lambda_i}{\lambda_1},i=2,3$ tienden a $0$ y por lo tanto $v{(k)}$ tiende a $\frac{w_1}{||w_1||_2}$.
>[!tip] Comportamiento de autovalores complejos en el método de la potencia
>Si $\lambda_1$ no es un número real positivo, $\frac{\lambda_1}{|\lambda_1|}$ es un número (complejo) de módulo 1, y $\left(\frac{\lambda_1}{|\lambda_1|}\right)^k$ tiene módulo 1 para todo $k$ pero no tiene límite. Sin embargo, a partir de esta sucesión podemos aproximar el valor de $\lambda_1$ y conociendo $\lambda_1$ podemos obtener fácilmente un autovector correspondiente a este autovalor. Esta estrategia para calcular el autovalor de módulo máximo se conoce como el **método de la potencia**. Vamos a presentar el método realizando algunas suposiciones que facilitan la demostración. El método se puede aplicar en condiciones más generales, incluso cuando $A$ no es diagonalizable, pero en ese caso necesitamos utilizar la forma de Jordan de la matriz u otras herramientas que escapan el contenido de este apunte para la demostración.
## Lema del cociente de Rayleigh
>[!note] Lema del cociente de Rayleigh para autovectores
>Si $v$ es un autovector de $A\in\mathbb{C}^{n\times n}$ de autovalor $\lambda$, el **cociente de Rayleigh** $$r(v)=\frac{v^TAv}{v^Tv}$$es igual al autovalor $\lambda$.
### Demostración
La demostración es inmediata usando que si $v$ es autovector de autovalor $\lambda$, $Av=\lambda v$. A partir de la definición tenemos muchas formas de calcular el autovalor correspondiente a un autovector. Se puede demostrar que cuando tenemos una aproximación $v$ de un autovector, $r(v)$ es el valor que más se parece a un autovalor, en el sentido de que $\alpha=r(v)$ minimiza la norma de $||Av-\alpha v||_2$.
## Algoritmo 1: Método de la potencia
$A\in\mathbb{C}^{n\times n},v^{(0)}\in\mathbb{C}^n;$
$\text{for }k=1,2,...\ \text{do}$
	$v^{(k)}=\frac{Av^{(k-1)}}{||Av^{(k-1)}||_2}$
	$r_k)\frac{(v^{(k)})^TAv^{(k)}}{(v^{(k)})^Tv^{(k)}}$
$\text{end}$
## Convergencia Asintótica
>[!note] Teorema de Convergencia Asintótica
>Sea $A\in\mathbb{C}^{n\times n}$ una matriz diagonalizable con base de autovectores $\mathcal{B}=\{w_1,...,w_n\}$ y autovalores correspondientes $\{\lambda_1,...,\lambda_n\}$ con $|\lambda_1|>|\lambda_2|\geq\cdots\geq|\lambda_n|\geq0$ (existe un único autovalor de módulo máximo). Dado un vector $v^{(0)}\in\mathbb{C}^n$ tal que en su escritura en la base $\mathcal{B}$ el coeficiente de $w_1$ es no nulo, el algoritmo 1 satisface $|r_k-\lambda_1|\to0$ cuando $k\to\infty$.
### Demostración
Por hipótesis, podemos escribir
$$v^{(0)}=a_1w_1+\cdots+a_nw_n$$
con $a_1\neq0$. Siguiendo los pasos del ejemplo, obtenemos
$$v^{(k)}=\left(\frac{\lambda_1}{|\lambda_1|}\right)^k\frac{a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+\cdots+a_n\left(\frac{\lambda_n}{\lambda_1}\right)^kw_n}{\left|\left|a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+\cdots+a_n\left(\frac{\lambda_n}{\lambda_1}\right)^kw_n\right|\right|_2}$$
Despejando,
$$\left(\frac{|\lambda_1|}{\lambda_1}\right)^kv^{(k)}=\frac{a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+\cdots+a_n\left(\frac{\lambda_n}{\lambda_1}\right)^kw_n}{\left|\left|a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+\cdots+a_n\left(\frac{\lambda_n}{\lambda_1}\right)^kw_n\right|\right|_2}$$
y por lo tanto, llamando $\tilde{v}^{(k)}=\left(\frac{|\lambda_1|}{\lambda_1}\right)^kv^{(k)}$,
$$\tilde{v}^{(k)}\stackrel{k\to\infty}{\longrightarrow}\frac{a_1w_1}{||a_1w_1||_2}$$
Como $\tilde{v}^{(k)}$ es igual a $v^{(k)}$ multiplicando por un escalar, podemos fácilmente ver que
$$r(\tilde{v}^{(k)})=\frac{(\tilde{v}^{(k)})^TA\tilde{v}^{(k)}}{(v^{(k)})^Tv^{(k)}}=\frac{(v^{(k)})^TAv^{(k)}}{(v^{(k)})^Tv^{(k)}}=r(v^{(k)})$$
Como $\tilde{v}^{(k)}$ tiende a un autovector de $A$, concluimos que
$$\lim_{k\to\infty}r(v^{(k)})=\lim_{k\to\infty}r(\tilde{v}^{(k)})=\lambda$$
### Relación con el lema de Rayleigh y el método de la potencia
La conexión entre este teorema y el lema de Rayleigh se basa en la *cota* que el cociente de Rayleigh establece para los autovalores en el proceso iterativo.
El *lema de Rayleigh*, o cociente de Rayleigh, establece que dado un vector no nulo $v \in \mathbb{C}^n$ y una matriz $A$, el valor $r(v) = \frac{v^* A v}{v^* v}$ es una aproximación al autovalor asociado al autovector más cercano a $v$. En el contexto del método de la potencia, este cociente se usa iterativamente para aproximar el autovalor dominante $\lambda_1$ de una matriz $A$, aprovechando el hecho de que, al aplicar sucesivas potencias de $A$ sobre un vector inicial $v^{(0)}$ (que tiene una componente no nula en la dirección del autovector dominante $w_1$), las contribuciones de los autovectores asociados a los autovalores de menor módulo se atenúan.
Con cada iteración del método de la potencia, el cociente de Rayleigh calculado para el vector iterado $v^{(k)} = A^k v^{(0)}$ tiende a aproximarse a $\lambda_1$, debido a que la componente en la dirección del autovector dominante $w_1$ se amplifica exponencialmente en comparación con las componentes en otras direcciones. Esto explica por qué $|r_k - \lambda_1| \to 0$ cuando $k \to \infty$: el cociente de Rayleigh converge asintóticamente al autovalor dominante, proporcionando una cota de error en cada paso de la iteración.
En resumen, el *lema de Rayleigh* permite cuantificar la aproximación del cociente de Rayleigh al autovalor dominante durante el proceso iterativo del método de la potencia, justificando así la convergencia hacia $\lambda_1$ y proporcionando una base teórica para la precisión del método.
# Norma 2 de matrices
Dada una matriz $A\in\mathbb{K}^{n\times n}$ y una norma vectorial $||\cdot||_n$ en $\mathbb{K}^n$, definimos la norma inducida
$$||A||=\max_{0\neq x\in\mathbb{K}^m}\frac{||Ax||}{||x||}=\max_{x\in\mathbb{K}^m,||x||=1}||Ax||$$
## Norma 2 de matrices diagonales
En el caso de matrices diagonales, el cálculo de normas matriciales es muy simple.
>[!note] Norma 2 de matrices diagonales
>Si $A$ es diagonal $$\begin{array}{}A=\begin{pmatrix}\lambda_1 & 0 & 0 & \cdots & 0 \\ 0 & \lambda_2 & 0 & \cdots & 0 \\ 0 & 0 & \lambda_3 & \cdots & 0 \\ \vdots & \vdots & & \ddots & \vdots \\ 0 & 0 & \cdots & & \lambda_n\end{pmatrix}\\||A||_2=\max\{|\lambda_i|\}\end{array}$$
### Radio espectral
>[!note] Radio espectral
>Dada $A\in\mathbb{K}^{n\times n}$ definimos el radio espectral de $A$ como $$\rho(A)=\max\{|\lambda|,\ \text{con $\lambda$ autovalor de $A$}\}$$
### Norma 2 de matrices simétricas
Supongamos que $A\in\mathbb{R}^{n\times n}$ es simétrica, sea $\mathcal{B}=\{v_1,...,v_n\}$ una base ortonormal de autovectores y $\{\lambda_1,\lambda_2,\cdots,\lambda_n\}\subset\mathbb{R}$ los respectivos autovalores. Sea $x\in\mathbb{R}^n$ tal que $||x||_2=1$, escribiendo $x$ en la base $\mathcal{B}$ se tiene
$$x=\sum_{i=1}^n\alpha_i v_i$$
y por ser $\mathcal{B}=\{v_1,...,v_n\}$ un conjunto ortonormal tenemos que, considerando la siguiente observación:
>[!note] Norma de bases ortonormales
>Para toda base ortonormal $\mathcal{B}=\{v_1,...,v_n\}$, si $$w=\sum_{i=1}^n\alpha_iv_i$$entonces
>$$||w||^2=\sum_{i=1}^n|\alpha_i|^2,$$

tenemos que:
$$1=||x||_2^2=\sum_{i=1}^n\alpha_i^2$$
Por otro lado, como $Ax=\sum_{i=1}^n\alpha_i\lambda_iv_1$, usando nuevamente esta observación, tenemos
$$||Ax||_2^2=\sum_{i=1}^n\alpha_i^2\lambda_i^2$$
Por lo tanto,
$$\frac{||Ax||^2}{||x||^2}=\sum_{i=1}^n\alpha_i^2\lambda_i^2\leq\rho(A)^2,\qquad\forall x,||x||\leq1$$
ya que cada $|\lambda_i|\leq\rho(A)$. Por otro lado, llamando $j$ a un índice tal que $|\lambda_j|=\rho(A)$ y considerando el vector de norma 1 $x=v_j$, se tiene
$$||Ax||^2=|\lambda_j|^2||x||_2^2=\rho(A)^2$$
En definitiva
$$\sup_{||x||=1}||Ax||^2=\rho(A)$$
Resultado que recuadramos debajo
>[!note] Norma 2 de una matriz simétrica
>Si $A\in\mathbb{R}^{n\times n},A^T=A$ entonces
>$$||A||_2=\rho(A)$$
## Norma 2 generalizada
Sea $A\in\mathbb{R}^{n\times n}$ no necesariamente simétrica. Queremos calcular
$$\sup_{||x||_2=1}||Ax||_2^2=(Ax)^T=x^T(A^TA)x$$
Si llamamos $M=A^TA$, resulta simétrica. Consideremos una base $\mathcal{B}=\{v_1,...,v_n\}$ ortonormal de autovectores de $M$ y $\{\lambda_1,\lambda_2,...,\lambda_n\}\subset\mathbb{R}$ los respectivos autovalores. Notemos que $M$ es semidefinida positiva, pues la ecuación previa se ve que para todo $x\in\mathbb{R}^n$
$$x^TMx\geq0$$
y en particular se sigue que $\lambda_i\geq0$ pues si $v_i\in\text{ker}(A)$ debe ser $\lambda_i=0$ y en caso contrario
$$0<v_i^TMv_i=v_i^T\lambda_iv_i=\lambda_i$$
Con esta información escribimos como antes
$$x=\sum_{i=1}^n\alpha_iv_i,\quad1=||x||_2^2=\sum_{i=1}^n\alpha_i^2,\quad Mx=\sum_{i=1}^n\alpha_iMv_i=\sum_{i=1}^n\lambda_i\alpha_iv_i$$
de donde por la ortonormalidad de los $v_i$
$$x^TMx=\sum_{i=1}^n\lambda_i\alpha_i^2$$
Luego,
$$\sup_{||x||_2=1}||Ax||_2^2=\sup_{\sum_{i=1}^n\alpha_i^2=1}\sum_{i=1}^n\lambda_i\alpha_i^2=\lambda_\max=\rho(M)$$
En definitiva tenemos
>[!note] Norma 2
>Sea $A\in\mathbb{R}^{n\times n}$, entonces
>$$||A||_2=\sqrt{\rho(A^TA)}$$
### Ínfimo de las normas subordinadas
>[!note] Ínfimo de las normas subordinadas
>Sea $A\in\mathbb{K}^{n\times n}$, entonces
>$$\rho(A)=\inf\{||A||:||\cdot||\ \text{es una norma subordinada}\}$$
#### Demostración
Solo vamos a probar que $\rho(A)\leq||A||$ para cualquier norma subordinada. La otra desigualdad es más complicada y no la demostraremos aquí.
Si $\mathbb{K}=\mathbb{C}$ la demostración es más sencilla. Sea $\lambda$ autovalor de $A$ y sea $||\cdot||$ una norma matricial en $\mathbb{C}^{n\times n}$ la subordinada a la norma vectorial $||\cdot||$. Luego existe $0\neq v$ autovector de autovalor $\lambda$. Podemos suponer que $||v||=1$ y entonces
$$Av=\lambda v$$
de donde $||Av||\leq|\lambda|\,||v||=|\lambda|$. Por lo tanto, para todo $\lambda$ autovalor
$$||A||=\sup_{||x||=1}||Ax||\geq||Av||\geq|\lambda|$$
de donde
$$||A||\geq\rho(A)$$
Si $\mathbb{K}=\mathbb{R}$ hay un problema al repetir el argumento previo, donde vemos que funciona si $\lambda\in\mathbb{R}$, sin embargo el autovalor y el autovector podrían ser complejo y no tener sentido la expresión $||v||$ si se trata de una norma en $\mathbb{R}^n$. Supongamos entonces que $\lambda=a+ib\in\mathbb{C}$ y escribimos entonces,
$$v=\text{Re}(v)+i\,\text{Im}(v)$$
y
$$Av=(a+ib)(\text{Re}(v)+i\,\text{Im}(v))$$
luego
$$A\,\text{Re}(v)=a\,\text{Re}(v)-b\,\text{Im}(v)\qquad A\,\text{Im}(v)=b\,\text{Re}(v)+b\,\text{Im}(v)$$
y entonces
$$\begin{array}{}
||A\,\text{Re}(v)||_2^2=a^2||\text{Re}(v)||_2^2-2ab\left<\text{Im}(v),\text{Re}(v)\right>+b^2||\text{Im}(v)||_2^2\\
||A\,\text{Re}(v)||_2^2=b^2||\text{Re}(v)||_2^2-2ab\left<\text{Im}(v),\text{Re}(v)\right>+a^2||\text{Im}(v)||_2^2
\end{array}$$
sumando miembro a miembro
$$||A\,\text{Re}(v)||_2^2+||A\,\text{Im}(v)||_2^2=(a^2+b^2)\left(||\text{Re}(v)||_2^2+||\text{Im}(v)||_2^2\right)$$
de donde
$$\left(||\text{Re}(v)||2^2+||\text{Im}(v)||_2^2\right)||A||_2^2\geq(a^2+b^2)\left(||\text{Re}(v)||_2^2+||\text{Im}(v)||_2^2\right)$$
y entonces
$$||A||_2^2\geq(a^2+b^2)=|\lambda|^2$$
Es decir -ya que $\lambda$ es arbitriario-
$$||A||_2\geq\rho(A)$$
Sea ahora $||\cdot||$ una norma arbitraria de $\mathbb{R}^n$ y $k\in\mathbb{N}$: por un lado sabemos que $\rho(A)^k\leq\rho(A^k)$ y que además (por la equivalencia de normas), existe una constante $C$-independiente de $k$- tal que $||A^k||_2\leq C||A^k||$. Por lo tanto,
$$\rho(A)^k\leq\rho(A^k)\leq||A^k||_2\leq C||A^k||\leq C||A||^k$$
donde hemos usado la desigualdad de normas subordinadas en la última desigualdad. Tomando raíz $k$-ésima
$$\rho(A)\leq C^{\frac{1}{k}}||A||$$
tomando $k\to\infty$ se obtiene $\rho(A)\leq||A||$, como queríamos probar.
### Condición de Convergencia a Cero de las Potencias de una Matriz
Como aplicación, obtenemos un criterio para determinar cuándo las potencias de una matriz tienden a la matriz nula.
>[!note] Condición de Convergencia a Cero de las Potencias de una Matriz
>Sea $M\in\mathbb{K}^{n\times n}$, se tiene que
>$$\lim_{k\to\infty}||M^k||=0$$
>para cierta norma $||\cdot||$ (y por ende para cualquier otra ya que son todas equivalentes) si y solo si $\rho(M)<1$.
#### Demostración
Si $\rho(M)<1$, entonces existe $||\cdot||$ tal que $||M||=\mu<1$ y listo.
Por otro lado,
$$\rho(M)^k=\rho(M^k)\leq||M^k||\to0$$
de donde necesariamente $\rho(M)<1$.
# Factorización QR
>[!abstract] Factorización QR
>Algoritmo que nos permite aproximar todos los autovalores y autovectores de una matriz de forma muy eficiente.
## Ideas principales
La demostración completa de este algoritmo excede lo visto, por lo que solo veremos ideas principales.
### Extensión del método de la potencia a subespacios
El método de la potencia nos permite calcular el autovalor de mayor módulo mediante una sucesión $\{v^{(k)}\}$ con
$$v^{(k)}=\left(\frac{\lambda_1}{|\lambda_1|}\right)^k\frac{a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+\cdots+a_n\left(\frac{\lambda_n}{\lambda_1}\right)^kw_n}{\left|\left|a_1w_1+a_2\left(\frac{\lambda_2}{\lambda_1}\right)^kw_2+\cdots+a_n\left(\frac{\lambda_n}{\lambda_1}\right)^kw_n\right|\right|_2}$$
Si bien, en general, $\{v^{(k)}\}\not\to w_1$ (ni a un múltiplo fijo), podemos intuitivamente pensar en convergencia de subespacios:
$$\left<v^{(k)}\right>\to\left<w_1\right>$$
Para calcular el resto de los autovalores y autovectores aplicamos el método (sin normalizar) a dos vectores en vez de uno.
$$\begin{array}{}
	v_1^{(k)}=(\lambda_2)^k\left(a_1\left(\frac{\lambda_1}{\lambda_2}\right)^kw_1+a_2w_2+a_3\left(\frac{\lambda_3}{\lambda_2}\right)^kw_3+\cdots+a_n\left(\frac{\lambda_n}{\lambda_2}\right)^kw_n\right)\\
	v_2^{(k)}=(\lambda_2)^k\left(b_1\left(\frac{\lambda_1}{\lambda_2}\right)^kw_1+b_2w_2+b_3\left(\frac{\lambda_3}{\lambda_2}\right)^kw_3+\cdots+b_n\left(\frac{\lambda_n}{\lambda_2}\right)^kw_n\right)
\end{array}$$
Intuitivamente vemos que $\left<v_1^{(k)},v_2^{(k)}\right>\to\left<w_1,w_2\right>$ cuando $k\to\infty$.
1. Si bien $\left<v_1^{(k)},v_2^{(k)}\right>\to\left<w_1,w_2\right>$, individualmente $\left<v_1^{(k)}\right>\to\left<w_1\right>$ y $\left< v_2^{(k)}\right>\to\left<w_2\right>$.
2. Si logramos "sacarle" a $v_2^{(k)}$ la componente en la dirección de $w_1$, vamos a obtener $\left<\tilde{v}_2^{(k)}\right>\to\left<w_2\right>$.
3. Si $\{w_1,..,w_n\}$ base ortonormal y $v=a_1w_1+\cdots+a_nw_n$, tenemos $(v,w_1)=a_1$.
4. Por lo tanto, si conocemos $w_1$, podemos hacer $v-(v,w_1)w_1=a_2w_2+\cdots+a_nw_n$.
5. En nuestro caso, $\left<v_2^{(k)}-(v_2^{(k)},w_1)w_1\right>\to\left<w_2\right>$.
En general, no conocemos $w_1$ pero podemos usar $v_1^{(k)}$.
Si $\left<v_1^{(k)}\right>\to\left<w_1\right>$ y definimos $u_1^{(k)}=v_1^{(k)},\tilde{u}_1^{(k)}=\frac{u_1^{(k)}}{||u_1^{(k)}||}$,
$$u_2^{(k)}=v_2^{(k)}-\left<\tilde{u}_1^{(k)},v_2^{(k)}\right>\tilde{u}_1^{(k)}\quad\text{y}\quad\tilde{u}_2^{(k)}=\frac{u_2^{(k)}}{||u_2^{(k)}||}$$
obtenemos $\left<\tilde{u}_2^{(k)}\right>\to\left<w_2\right>$.
#### Observación
$\{\tilde{u}_1^{(k)},\tilde{u}_2^{(k)}\}$ es el resultado de aplicar [[Teórica 3 - Sistemas lineales y factorización#Método de ortonormalización de Gram-Schmidt|Gram-Schmidt]] a $\{v_1^{(k)},v_2^{(k)}\}$.
### Aplicación en matrices simétricas
Si extendemos el razonamiento a un conjunto inicial $\{v_1^{(0)},...,v_n^{(0)}\}$ de vectores linealmente independientes podemos generar sucesiones de vectores tales que $\left<\tilde{u}_i^{(k)}\right>\to\left<w_i\right>$.
Trabajando con matrices, podemos aplicar el método simultáneamente a un conjunto de vectores. Obtenemos el siguiente método para matrices simétricas.

>[!info] Algoritmo 2: Iteración simultánea
>$A$ matriz (simétrica) con autovalores $|\lambda_1|>|\lambda_2|>...>|\lambda_n|$;
>$B^{(0)}=I_n$ (o cualquier matriz de rango máximo);
>$\text{for}\ k=1,2,...\ \text{do}$
>	$\tilde{B}^{(k)}=AB^{(k-1)}$;
>	$\tilde{B}^{(k)}=Q^{(k)}R^{(k)}$, descomposición $QR$ de $\tilde{B}^{(k)}$;
>	$\tilde{B}^{(k)}=Q^{(k)}$, el resultado de aplicar [[Teórica 3 - Sistemas lineales y factorización#Método de ortonormalización de Gram-Schmidt|Gram-Schmidt]] a las columnas de $\tilde{B}^{(k)}$
>$\text{end}$

Si $c_i^{(k)}$ es la $i$-ésima columna de $B^{(k)}$,
$$\left<c_i^{(k)}\right>\to\left<w_i\right>$$
Para aproximar los autovalores, utilizamos la matriz $A$ en la base $\mathcal{B}^{(k)}$ dada por $B^{(k)}$, y obtenemos
$$[A]_{\mathcal{B}^{(k)}}=(B^{(k)})^TAB^{(k)}$$
la cual es una matriz diagonal con los autovalores en la diagonal.
### Extensión a matrices generales
Para aplicar el método anterior en matrices generales, con autovalores $|\lambda_1|>|\lambda_2|>...>|\lambda_n|$, tenemos que considerar que, si bien no se cumple $\left<\tilde{v}_i^{(k)}\right>\to\left<w_i\right>$ para todo $i$, obtenemos en cambio
$$\begin{cases}
	\left<\tilde{v}_1^{(k)}\right> & \to& \left<w_1\right>\\
	\left<\tilde{v}_1^{(k)},\tilde{v}_2^{(k)}\right> & \to& \left<w_1,w_2\right>\\
	 & \vdots \\
	\left<\tilde{v}_1^{(k)},...,\tilde{v}_n^{(k)}\right> & \to & \left<w_1,...,w_n\right>
\end{cases}$$
Por lo tanto,
$$[A]_{\mathcal{B}^{(k)}}=(B^{(k)})^TAB^{(k)}$$
tiende a una matriz triangular superior con los autovalores en la diagonal.
Finalmente, acomodando las operaciones convenientemente, obtenemos el siguiente método
#### Algoritmo 3: Método QR
$A^{(0)}=A\in\mathbb{C}^{n\times n}$ con autovalores $|\lambda_1|>|\lambda_2|>...>|\lambda_n|$;
$\text{for}\ k=0,1,2,...\ \text{do}$
	Calcular la descomposición $QR$ de $A^{(k)},A^{(k)}=Q^{(k)}R^{(k)}$;
	$A^{(k+1)}=R^{(k)}Q^{(k)}$
$\text{end}$

>[!note] Convergencia de matrices a triangulares superiores mediante Factorización QR
>La matriz $A^{(k)}$ tiende a una matriz triangular superior con los autovalores de $A$ en la diagonal.

>[!tip] Invarianza de Autovalores en la Secuencia de Matrices mediante Factorización QR
>La matriz que se obtiene en cada paso es semejante a la matriz del paso anterior:
>$$(Q^{(k)})^TA^{(k)}Q^{(k)}=(Q^{(k)})^T(Q^{(k)}R^{(k)})Q^{(k)}=R^{(k)}Q^{(k)}=A^{(k+1)}$$
>Por lo tanto, las matrices $A^{(k)}$ tienen los mismos autovalores para todo $k=0,1,2,...$
# Teorema de los círculos de Gershgorin
Esta es otra herramienta en la caja de los buscadores de autovalores.
Dada $A\in\mathbb{C}^{n\times n}$, definimos para $1\leq i\leq n,R_i$ la suma de los módulos de los elementos de la $i$-ésima fila de $A$ fuera de la diagonal,
$$R_i=\sum_{j\neq i}|a_{ij}|$$
Llamamos $D(a_{ii},R_i)\subset\mathbb{C}$ al disco cerrado centrado en $a_{ii}$ de radio $R_i$.

>[!note] Localización de Autovalores mediante los Discos de Gershgorin 
>Todos los autovalores de $A$ están en alguno de los discos de Gershgorin $D(a_{ii},R_i)$.