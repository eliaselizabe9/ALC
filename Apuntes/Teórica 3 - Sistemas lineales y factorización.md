# Ortonormalización
## Método de ortonormalización de Gram-Schmidt
>[!note] Método de ortonormalización de Gram-Schmidt
>Sea $(V,<,>)$ un espacio vectorial con ==producto interno definido== y sea $\{v_1,...,v_n\}$ una base de $V$. Existe una base ortonormal $B=\{w_1,...,w_n\}$ de $V$ tal que
>$$<v_1,...,v_k>=<w_1,...,w_k>\qquad\forall 1\leq k\leq n$$
### Demostración
Se construirá una base ortogonal $\{z_1,...,z_n\}\in V$ que cumpla
$$<z_1,...,z_k>=<v_1,...,v_k>\qquad\forall 1\leq k\leq n$$
Normalizando los vectores de esta base se obtendrá la base ortonormal buscada.
Construiremos los vectores de la base recursivamente.
- Tomamos $z_1=v_1$, que satisface $<z_1>=<v_1>$. Buscamos $z_2\in V$ con $<z_2,z_1>=0$ y tal que $<z_1,z_2>=<v_1,v_2>$. Esta segunda condición vale si y sólo si $z_2$ es de la forma $z_2=a\cdot v_1+b\cdot v_2$ con $b\neq0$. Podemos suponer entonces que $b=1$, es decir, $z_2=v_2+a\cdot v_1$ y buscar $a$ de manera que se cumpla la primera condición: $$0=<z_2,z_1>=<v_2+av_1,v_1>=<v_2,v_1>+a<v_1,v_1>$$lo que implica que $a=\frac{-<v_2,v_1>}{||v_1||^2}$. Luego el vector $$z:2=v_2-\frac{<v_2,v_1>}{||v_1||^2}\cdot v_1=v_2-\frac{<v_2,z_1>}{||z_1||^2}\cdot z_1$$satisface las condiciones requeridas.
- Supongamos construidos $z_1,...,z_r\in V$ tales que
	- $<z_i,z_j>=0$ si $i\neq j$.
	- $<z_1,...,z_k>=<v_1,...,v_k>\quad\forall 1\leq k\leq r$
	Consideramos el vector $$z_{r+1}=v_{r+1}-\sum_{i=1}^r\frac{<v_{r+1},z_i}{||z_i||^2}\cdot z_i$$Se tiene que:
	- $<z_1,...,z_r,z_{r+1}=<z_1,...,z_r,v_{r+1}>=<v_1,...,v_r,v_{r+1}>$
	- para cada $j\neq r$: $$\begin{array}{lll}<z_{r+1},z_j> & = &\left<v_{r+1}-\sum_{i=1}^r\frac{<v_{r+1},z_i>}{||z_i||^2}\cdot z_i,z_j\right>\\ & = & <v_{r+1},z_j>-\sum_{i=1}^r\frac{<v_{r+1},z_i>}{||z_i||^2}\cdot <z_i,z_j>\\ & = & <v_{r+1},z_j>-\frac{<v_{r+1},z_j>}{||z_j||^2}\cdot <z_j,z_j>\\ & = & 0\end{array}$$Luego, $z_{r+1}$ satisface las condiciones requeridas.
De este modo, al concluir el $n$-ésimo paso se obtiene una base ortogonal $\{z_1,...,z_n\}$ de $V$ que además satisface $<z_1,...,z_k>=<v_1,...,v_k>$ para cada $1\leq k\leq n$.
Finalmente, para cada $1\leq i\leq n$ consideramos al vector $w_i=\frac{z_i}{||z_i||}$. Luego, el conjunto $B=\{w_1,...,w_n\}$ resulta una base ortonormal de $V$ que cumple lo pedido.
### Corolario de la Extensión de la Base Ortonormal en Subespacios
>[!note] Extensión de la Base Ortonormal en Subespacios
>Sea $(V,<,>)$ un espacio vectorial de dimensión finita con producto interno, y sea $S$ un subespacio de $V,S\neq\{0\}$. Entonces existe una base ortonormal de $V$ que contiene una base ortonormal de $S$.
#### Demostración
Sea $\{s_1,...,s_r\}$ una base de $S$. Existen $v_{r+1},...,v_n\in V$ tales que $B=\{s_1,...,s_r,v_{r+1},...,v_n\}$ es una base de $V$. Aplicando el proceso de Gram-Schmidt a $B$ se obtiene una base ortonormal $B'=\{w_1,...,w_r,w_{r+1},...,w_n\}$ de $V$ que satisface
$$<w_1,...,w_r>=<s_1,...,s_r>=S$$
En consecuencia, $\{w_1,...,w_r\}$ es una base ortonormal de $S$ incluida en la base ortonormal $B'$ de $V$.
### Ejemplo
Dada la base $B=\{(1,0,i),(1,1,2+i),(0,0,1)\}$ de $\mathbb{C}^3$, ortonormalizarla usando el método de Gram-Schmidt.
Notaremos $v_1=(1,0,i)$, $v_2=(1,1,2+i)$ y $v_3=(0,0,1)$.
- Comenzamos tomando $z_1$ como el primer vector de la base: $z_1=(1,0,i)$.
- Construimos ahora $z_2$ aplicando la fórmula $$z_{r+1}=v_{r+1}-\sum_{i=1}^r\frac{<v_{r+1},z_i>}{||z_i||^2}\cdot z_i$$para $r=1$: $$\begin{array}{lll}z_2 & = & v_2-\frac{<v_2,z_1>}{||z_1||^2}\cdot z_1=(1,1,2+i)-\frac{<(1,1,2+i),(1,0,i)>}{||(1,0,i)||^2}\cdot(1,0,i)\\ & = & (1,1,2+i)-(1-i)\cdot(1,0,i)=(i,1,1)\end{array}$$
- Finalmente, hallamos $z_3$ aplicando nuevamente la fórmula vista previamente para $r=2$:$$\begin{array}{lll}z_3 & = & v_3-\frac{<v_3,z_1>}{||z_1||^2}\cdot z_1-\frac{<v_3,z_2>}{||z_2||^2}\cdot z_2\\ & = & (0,0,1)-\frac{<(0,0,1),(1,0,i)>}{||(1,0,i)||^2}\cdot(1,0,i)-\frac{<(0,0,1),(i,1,1)>}{||(i,1,1)||^2}\cdot(i,1,1)\\ & = & (0,0,1)+\frac{i}{2}\cdot(1,0,i)-\frac{1}{3}\cdot(i,1,1)\\ & = & \left(\frac{i}{6},\frac{-1}{3},\frac{1}{6}\right)\end{array}$$
El conjunto $\{z_1,z_2,z_3\}$ es una base ortogonal de $\mathbb{C}^3$. Dividiendo cada vector por su norma, obtenemos
$$\begin{array}{lll}
	w_1 & = & \frac{z_1}{||z_1||}=\frac{1}{\sqrt{2}}\cdot(1,0,i)=\left(\frac{1}{\sqrt{2}},0,\frac{i}{\sqrt{2}}\right)\\ w_2 & = & \frac{z_2}{||z_2||}=\frac{1}{\sqrt{3}}\cdot(i,1,1)=\left(\frac{i}{\sqrt{3}},\frac{1}{\sqrt{3}},\frac{1}{\sqrt{3}}\right)\\ w_3 & = & \frac{z_3}{||z_3||}=\sqrt{6}\left(\frac{i}{6},\frac{-1}{3},\frac{1}{6}\right)=\left(\frac{\sqrt{6}i}{6},\frac{-\sqrt{6}}{3},\frac{\sqrt{6}}{6}\right)
\end{array}$$
tales que $\{w_1,w_2,w_3\}$ es una base ortonormal de $\mathbb{C}^3$ que cumple $<v_1>=<w_1>$ y $<v_1,v_2>=<w_1,w_2>$.
## Proyección sobre un subespacio normalizado
Otra forma equivalente de escribir la proyección ortogonal de $w$ sobre el subespacio generado por $v$ es
$$P_v(w)=\frac{v}{||v||_2}\frac{v^*}{||v||_2}w$$
La matriz resultante, asociada a la proyección $P_v$, tiene rango 1. Las columnas resultantes son todas múltiplos de $v$.
Cuando $v$ es de norma uno
$$P_v(w)=vv*w$$
Considerando esto y siendo que $w=\sum_{i=1}^n\alpha_iv_i$. Si en vez de una base completa solo contamos con un conjunto de vectores ortonormales $\mathcal{C}=\{v_1,\cdots,v_k\}\subset\mathbb{K}^n,k<n$, y llamamos $S=<v_1,\cdots,v_k>$ al subespacio generado por estos vectores, vemos que la expresión
$$P_S(w)=\sum_{1\leq i\leq k}P_{v_i}(w)$$
verifica
$$(w-P_S(w))\perp v_j,\qquad\forall 1\leq j\leq k$$
pues
$$v_j^*(w-P_S(w))=v_j^*w-\sum_{1\leq i\leq k}v_j^*P_{v_i}(w)=v_j^*w-v_j^*w=0$$
Luego, por ser ortogonal a todos los generadores de $S$ lo será a todos los elementos de $S$. En otras palabras, ==teniendo una base ortonormal de $S$ es muy fácil construir la proyección ortogonal de un elemento $w$ sobre $S$ usando la expresión anterior==. Como además, cada $P_{v_i}$ admite una representación matricial, lo mismo ocurre con $P_S$.
![[Pasted image 20240427235453.png]]
En efecto, tomando $w\in\mathbb{K}^n$ tenemos, llamando $M$ a la matriz anterior que
![[Pasted image 20240427235539.png]]
es decir que coincide con la expresión de proyección.
# Factorización LU
>[!note] Factorización LU
>Escribir una matriz invertible $A$ como producto de dos factores:
>- **Triangular inferior** $L$
>- **Triangular superior** $U$
>$$L=\begin{pmatrix}* & 0 & \cdots & 0 \\ * & * & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ * & * & \cdots & *\end{pmatrix}, U=\begin{pmatrix}* & * & \cdots & * \\ 0 & * & \cdots & * \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & *\end{pmatrix}$$
>y
>$$A=LU$$

La ventaja práctica de este sistema radica en que para resolverlo basta con resolver los sistemas
$$Ly=b\qquad\land\qquad Ux=y$$
y cada sistema triangular puede resolverse por sustitución (regresiva o progresiva) en $\sim n^2$ operaciones.
## Descomposición $LU=A$
### Lema de Inversión de Matriz Modificada
>[!note] Lema de Inversión de Matriz Modificada
>Sea $u,v\in\mathbb{R}^n$, tales que $v^tu\neq1$ entonces $I-uv^t$ es invertible y además
>$$(I-uv^t)^{-1}=I+\frac{uv^t}{1-v^tu}$$
#### Demostración
$$(I-uv^t)\left(I+\frac{uv^t}{1-v^tu}\right)=I+uv^t\left(\frac{1}{1-v^tu}-1\right)-\frac{uv^tuv^t}{1-v^tu}=I$$
#### Costo algorítmico de invertir $M_i$
Considerando el lema anterior y observando que $M_i=I-z_ie_i^t$ tenemos que $e_i^tz_i=0$, resultando en que
$$M_i^{-1}=I+z_ie_i^t$$
de donde
$$L=M_1^{-1}M_2^{-1}\cdots M_{(n-1)}^{-1}=(I+z_1e_1^t)(I+z_2e_2^t)\cdots(I+z_{n-1}e_{n-1}^t)$$

>[!info] Costo de invertir la matriz $M_i$
>El costo de invertir $M_i$ es casi nulo, apenas un **cambio de signos debajo de la diagonal**.

Además de la sencillez computacional de $M_i^{-1}$, podemos expresar $L$ sin cálculos adicionales. Como $e_i^tz_k=0$ no solo para $k=i$ sino para todo $i\leq k\leq n$, se puede desarrollar la expresión para $L$ y obtener
$$L=I+\sum_{i=1}^{n-1}z_ie_i^t$$
o dicho de otra forma
$$L=\begin{pmatrix} 1 & 0 & \cdots & 0 & 0 \\ \frac{a_{2,1}}{a_{1,1}} & 1 & \cdots & 0 & 0 \\ \frac{a_{3,1}}{a_{1,1}} & \frac{a_{3,2}^{(1)}}{a_{2,2}} & \vdots & \ddots & \vdots \\ \vdots & \vdots & \ddots & 1 & \vdots \\ \frac{a_{n,1}}{a_{1,1}} & \frac{a_{n,2}^{(1)}}{a_{2,2}^{(1)}} & \cdots & \frac{a_{n,n-1}^{(n-1)}}{a_{n-1,n-1}^{(n-1)}} & 1\end{pmatrix}$$
o sea, **basta con almacenar los multiplicadores en sus respectivos lugares cambiando los signos**.
Una condición básica para la descomposición $LU$ es que todos los pivots sean no nulos. Eso no puede garantizarse bajo la única hipótesis de que $\text{det}(A)\neq0$ como se ve en el siguiente ejemplo.
### Criterio 1 de existencia de la factorización LU
Sea $A=\begin{pmatrix}0 & 1 \\ 1 & 0\end{pmatrix}$. Sin embargo, puede garantizarse si se requiere que todos los menores de la matriz $A$, sean invertibles.
>[!note] Condición de Invertibilidad para la existencia de la Factorización LU
>Sea $A\in\mathbb{K}^{n\times n}$, entonces existe factorización $LU$ de $A$ con $L,U$ invertibles sí y solo sí para todo $1\leq k\leq n,\text{det}(1:k,1:k))\neq0$ (i.e. **todos los determinantes de las submatrices de $A$ deben ser no nulos**).
#### Demostración
##### $\Longleftarrow$
Como $a_{11}=A(1,1)$ resulta $a_{11}\neq0$ y el primer pivot es no nulo. Realizada la eliminación gaussiana, vemos que se puede escribir
$$A=M_1^{-1}U_1$$
donde
$$M_1^{-1}U_1=
\begin{pmatrix}
	1 & 0 & \cdots & 0 & 0 \\
	\frac{a_{2,1}}{a_{1,1}} & 1 & \cdots & 0 & 0 \\
	\frac{a_{3,1}}{a_{1,1}} & 0 & \vdots & \ddots & \vdots \\
	\vdots & \vdots & \ddots & 1 & \vdots \\
	\frac{a_{n,1}}{a_{1,1}} & 0 & \cdots & 0 & 1
\end{pmatrix}
\begin{pmatrix}
	a_{1,1} & a_{1,2} & a_{1,3} & \cdots & a_{1,n} \\
	0 & a_{2,2}^{(1)} & a_{2,3}^{(1)} & \cdots & a_{2,n}^{(1)} \\
	0 & 0 \\
	\vdots & \vdots & & A^{(2)} & \\
	0 & 0 
\end{pmatrix}$$
Por lo tanto
$$A(1:2,1:2)=\begin{pmatrix}1 & 0 \\ \frac{a_{2,1}}{a_{1,1}} & 1\end{pmatrix}\begin{pmatrix}a_{1,1} & a_{1,2} \\ 0 & a_{2,2}^{(1)}\end{pmatrix}$$
y como por hipótesis $0\neq\text{det}(A(1:2,1:2))$ tenemos $\text{det}(A(1:2,1:2))=a_{1,1}a_{2,2}^{(1)}\neq0$, por lo tanto el segundo pivot no se anula. La demostración se sigue por inducción.
##### $\implies$
Notamos que si $A=LU$ con $L$ y $U$ invertibles, deben ser $l_{kk}\neq0\neq u_{kk}$ para todo $1\leq k\leq n$ ya que son triangulares. En particular, $\text{det}(L(1:k,1:k))\neq0\neq\text{det}(1:k,1:k))$ y por ser una de ellas triangular inferior y la otra superior, resulta $A(1:k,1:k)=L(1:k,1:k)U(1:k,1:k)$ de donde $\text{det}(A(1:k,1:k))\neq0$.
### Criterio 2 de existencia de la factorización LU
>[!note] Condición de Diagonal Dominante para la existencia de factorización $LU$
>Sea $A\in\mathbb{R}^{n\times n},\mathbb{C}^{n\times n}\ \text{EDD}$ entonces $A$ es invertible.
#### Diagonal dominante
>[!note] Diagonal dominante
>Una matriz $A\in\mathbb{R}^{n\times n},\mathbb{C}^{n\times n}$ se dice **diagonal dominante** (estrictamente diagonal dominante) y se denota $\text{DD}(\text{EDD})$ si y sólo si para todo $i,1\leq i\leq n$. $$\sum_{1\leq j\leq n,j\neq i}|a_{i,j}|\leq(z)|a_{i,i}|$$
#### Demostración
Supongamos que no. Existe entonces $0\neq v\in\mathbb{C}^n$ tal que $Av=0$. Sea $i,1\leq i\leq n$ tal que $0\neq|v_i|=||v||_{\infty}$, entonces
$$\sum_{j=1}^na_{i,j}v_j=0$$
de donde
$$|a_{i,i}|\leq\sum_{j=1,j\neq i}^n|a_{i,j}|\frac{|v_j|}{|v_i|}\leq\sum_{j=1,j\neq i}^n|a_{i,j}|$$
lo que contradice la $\text{EDD}$.
### Propiedad 3
>[!note] Propiedad de No-Ceros Pivote para Matrices $\text{EDD}$
>Si $A$ es $EDD$ entonces, trabajando con aritmética exacta, el proceso de eliminación de Gauss no produce pivots nulos.
#### Demostración
#### Observación
Las matrices $\text{DD}$ no son un artificio teórico, aparecen frecuentemente en la práctica como por ejemplo en la discretización de ecuaciones diferenciales.
### Teorema de Estabilidad de la Factorización LU
>[!note] Teorema de Estabilidad de la Factorización LU
>Sea $A\in\mathbb{K}^{n\times n}$, si el algoritmo $LU$ no produce pivots nulos, la factorización $LU$ producida por al máquina verifica
>$$\tilde{L}\tilde{U}=A+\Delta A$$
>donde $|\Delta A|\lesssim2(n-1)\epsilon(|A|+|\tilde{L}||\tilde{U}|)$
#### Observación
>[!tip] Controlar el tamaño de los factores $\tilde{L},\tilde{U}$
>Si el producto $|\tilde{L}||\tilde{U}|$ es muy grande no tendremos control (al menos teórico) sobre el error. El *pivoteo parcial* permite garantizar al menos el control $|L|\leq1$ a un costo razonable. Si bien, esto no necesariamente implica que $|\tilde{L}||\tilde{U}|$ sea pequeño, resulta suficiente en la mayoría de los casos prácticos.
## Descomposición $LU=PA$ (pivoteo parcial)
### Consideraciones iniciales de la eliminación gaussiana
Una matriz de permutaciones $P\in\mathbb{K}^{n\times n}$ es una matriz identidad por las filas permutadas. Alternativamente puede verse como una matriz cuyas filas están conformadas por los elementos de la base canónica no necesariamente ordenados. Es decir,
$$P=\begin{pmatrix}e_{i_1}^T\\----\\ e_{i_2}^T\\\vdots\\----\\e_{i_n}^T\end{pmatrix}$$
donde $i_1,i_2,...,i_n\in\{1,2,...,n\}$ son todos diferentes. Alternativamente podemos pensar en una matriz de permutaciones como la identidad con las columnas permutadas.
$$P=\begin{pmatrix}e_{i_1} & | & e_{i_2} & | & \cdots & | &e_{i_n}\end{pmatrix}$$
Notemos que si $v\in\mathbb{K}^n$, entonces
$$Pv=\begin{pmatrix}v_{i_1} \\ v_{i_2} \\ \vdots \\ v_{i_n}\end{pmatrix}\qquad\begin{pmatrix}v^T P = (v_{i_1},v_{i_2},\cdots,v_{i_n}\end{pmatrix}$$
Resumamos las siguientes propiedades de las permutaciones.
- Dada $A\in\mathbb{K}^{n\times n},P A$ coincide con $A$ pero con las filas permutadas y $AP$ coincide con $A$ pero con las columnas permutadas.
- Dada $x\in\mathbb{K}^n,Px$ coincide con $x$ pero con los elementos permutados.
- $\text{det}(P)=\pm1$
- $P^{-1}=P^T$
- Si $P_1$ y $P_2$ son permutaciones entonces $P=P_1P_2$ también lo es.
Si durante la eliminación de Gauss hallamos un pivot nulo, podemos intercambiar filas para continuar con el algoritmo (siempre que haya algún elemento no nulo con el cual proseguir). Mas aún, aunque el correspondiente pivot sea no nulo, podemos, casi al mismo costo, tomar el pivot más grande (en valor absoluto). Esa operación la podemos representar matricialmente multiplicando por una adecuada permutación.
### Procedimiento de pivoteo parcial
Seleccionamos el elemento con mayor valor absoluto en la primera columna de $A$. Dicho $a_{k,1}$, si el $\text{det}(A)\neq0$ tenemos que $a_{k,1}\neq0$. Permutamos la fila $k$ con la fila 1, lo cual se hace con una matriz $P_1$
$$P_1A=
\begin{pmatrix}
	a_{k,1} & a_{k,2} & \cdots & a_{k,n}\\
	a_{2,1} & a_{2,2} & \cdots & a_{2,n}\\
	\vdots & \vdots & \vdots & \vdots\\
	a_{1,1} & a_{1,2} & \cdots & a_{1,n}\\
	\vdots & \vdots & \vdots & \vdots\\
	a_{n,1} & a_{n,2} & \cdots & a_{n,n}
\end{pmatrix}$$
a partir de aquí construimos el multiplicador, $M_1$ de modo tal que
$$M_1P_1A=\begin{pmatrix}
	a_{k,1} & a_{k,2} & \cdots & a_{k,n} \\
	0 \\
	\vdots & & A^{(1)}\\
	0
\end{pmatrix}$$
Llegado a este punto, repetimos el procedimiento con la primer columna de $A^{(1)}$ que eventualmente requiere otra matriz $P_2$ antes de la construcción del nuevo multiplicador $M_2$. En definitiva,
$$M_{n-1}P_{n-1}\cdots M_1P_1A=U$$
se obtiene una nueva matriz triangular superior $U$, y donde $M_i=I-z_ie_i^T$ con
$$z_i=
\begin{pmatrix}
	0 \\ \vdots \\ 0 \\ z_{i+1}^{(i)} \\ \vdots \\ z_n^{(i)}
\end{pmatrix}$$
$|z_i|\leq1$ (o sea, los multiplicadores son menores o iguales a 1 en valor absoluto). El problema ahora es que en la serie de multiplicadores $M_nP_n$ que constituyen $U$ no parece verse el modo de reconstruir los factores $L$ y $P$ de manera sencilla. Sin embargo un hecho inesperadamente favorable resuelve esta cuestión. Notemos que por el orden de las iteraciones en el paso $k$ solo permutaremos en busca del mayor pivot (en caso de ser necesario), filas desde $k$ a $n$ (las filas previas no puesto que ya han sido procesadas). Entonces, $e_j^TP_k=e_j^T$ para todo $j<k$. En particular, si bien no es cierto que $P_k$ conmuta con $M_j$ sí es cierto, para todo $j<k$, que
$$P_kM_j=P_k(I-z_je_j^T)=(I-P_kz_je_j^T)P_k=\tilde{M}_jP_k$$
es decir que podemos pasar la permutación a la derecha permutando adecuadamente los multiplicadores de la matriz $M_j$. Lo importante es que esta nueva matriz, denotada con $\tilde{M}_j$ tiene la misma estructura que $M_j$. Así que volviendo a la sucesión de multiplicadores tenemos
$$M_{n-1}P_{n-1}\cdots M_1P_1A=\tilde{M}_{n-1}\cdots\tilde{M}_1P_{n-1}\cdots P_1 A=U$$
Argumentando como en el caso son pivoteo (usando que la estructura de las $\tilde{M}_i$ es igual a la de las $M_i$) y llamando $P=P_{n-1}\cdots P_1$ se llega a $PA=LU$.
### Propiedad 4
>[!note] Propiedad de pivots no nulos a partir del determinante
>Trabajando con aritmética exacta se tiene que si $\text{det}(A)\neq0$ el algoritmo de pivoteo parcial no tiene pivots nulos.
#### Observación
### Teorema de cota superior de la Norma de $U$
Controlar el tamaño de ambos factores $L$ y $U$ es posible. La idea se denomina **pivoteo completo**, esto implica hacer permutaciones de filas y de columnas para tomar siempre el mayor pivot posible. Así, en el paso inicial se busca el elemento $a_{k,l}$ con máximo valor absoluto y se permutan la fila y columna 1 con la fila y columna $k$ y $l$ respectivamente. En términos matriciales eso requiere de dos matrices de permutaciones $P_1,\tilde{P}_1$ antes de aplicar el paso de eliminación. Es decir
$$M_1,P_1A\tilde{P}_1=\begin{pmatrix}* & * & \cdots & * \\ = \\ \vdots & & A^{(1)}\\ 0\end{pmatrix}$$
Repitiendo el procedimiento se obtiene $U$
$$M_{n-1}P_{n-1}\cdots M_1P_1A\tilde{P}_1\cdots\tilde{P}_{n-1}=U$$
y argumentando como antes se obtienen ahora dos matrices de permutaciones $P,\tilde{P}$ tales que
$$PA\tilde{P}=LU$$
El pivoteo parcial no se lleva a cabo en general debido a su alto costo (calcule el número de comparaciones que debe realizar en cada) ya que hay otros métodos estables más económicos. En esta línea, Wilkinson probó el siguiente teorema.
>[!note] Teorema de cota superior para la Norma de $U$
>En aritmética exacta, para el algoritmo de pivoteo completo se tiene
>$$||U||_\infty\leq\sqrt{n}(2\cdot3^{1/2}\cdot4^{1/3}\cdots n^{1/(n-1)})^{1/2}||A||_\infty$$
# Factorización de Cholesky
En muchos casos es posible explotar la estructura de las matrices a la hora de resolver sistemas, con la finalidad de ahorrar cómputos o mejorar la estabilidad de los algoritmos.
## Teorema de Invertibilidad de matrices Definidas Positivas
>[!info] Teorema de Invertibilidad de Matrices Definidas Positivas
>Si $A$ es definida positiva entonces es invertible y también lo son todos sus menores.
### Demostración
En efecto, de no ser invertible, existiría $0\neq v$ tal que $Av=0$ lo que contradice la definida positividad ya que para ese $v$ se tendría $v*Av=0$.
Respecto de los menores el resultado se sigue del anterior. En efecto, sea $k$ fijo, $1\leq k\leq n$ y tomemos vectores $v\in\mathbb{K}^n$ de la forma $(\tilde{v},0,\cdots,0),\tilde{v}\in\mathbb{K}^k$ arbitrario. Resulta
$$0<v*Av=\tilde{v}*A(1:k,1:k)\tilde{v}$$
y de ahí la definida positividad de $A(1:k,1:k)$ y por ende su invertibilidad.
## Descomposición de Cholesky
### Existencia de la factorización de Cholesky
>[!note] Existencia de la factorización de Cholesky
>1. **Verificar que la matriz es definida positiva**: Si $A$ es definida positiva, entonces admite una descomposición $A=LU$.
>2. **Descomposición $LU$**: Podemos escribir $U=D\tilde{L}^*$, con $\tilde{L}$ triangular inferior con 1 en la diagonal y $D$ una matriz diagonal.
>3. **Obtener la matriz $L$**: De la ecuación $D(L^{-1}\tilde{L})^*=L^{-1}\tilde{L}D^*$, obtenemos que $L=\tilde{L}$.
>4. **Obtener la matriz $D$**: De la ecuación $D(L^{-1}\tilde{L})^*=L^{-1}\tilde{L}D^*$, obtenemos que $D=D^*$.
>5. **Factorización $LDL^*$**: En este punto, obtenemos la factorización $A=LDL^*$.
>6. **Verificar la definida positividad de $D$**: Como $A$ es definida positiva, entonces $D>0$.
>7. **Factorización de Cholesky**: Finalmente, obtenemos la factorización de Cholesky $A=LD^{\frac{1}{2}}D^{\frac{1}{2}}L^*=CC^*$, de donde $C$ es triangular con elementos diagonales positivos.
#### Explicación
Si $A$ es definida positiva admite una descomposición $A=LU$ (gracias a [[Teórica 3 - Sistemas lineales y factorización#Teorema de Invertibilidad de matrices Definidas Positivas|Teorema de Invertibilidad de matrices Definidas Positivas]] y [[Teórica 3 - Sistemas lineales y factorización#Diagonal dominante|Teorema de Diagonal Dominante]]) podemos escribir $U=D\tilde{L}^*$, con $\tilde{L}$ triangular inferior con 1 en la diagonal y $D$ una matriz diagonal. De ahí resulta
$$A=LD\tilde{L}=\tilde{L}D^*L^*$$
entonces
$$D(L^{-1}\tilde{L})^*=L^{-1}\tilde{L}D^*$$
y como el lado izquierdo de esta ecuación es triangular superior y el lado derecho es triangular inferior deberán ser ambos diagonales. Luego, debe serlo también $L{-1}\tilde{L}$. Por construcción, $\tilde{L}_{i,i}=1=L_{i,i}$, para todo $1\leq i\leq n$. Debido a la primera igualdad se tiene en particular que $L_{i,i}^{-1}$. De aquí resulta $L^{-1}\tilde{L}=I$, es decir $L=\tilde{L}$ y además $D=D^*$. En particular, se obtiene la factorización
$$A=LDL^*$$
Como además $A$ es definida positiva, entonces $D>0$, pues para todo $v\neq0$
$$0<v^*Av=v^*LDL^*v=(L^*v)^*DL^*v$$
de donde, llamando $w=L^*v$ y usando que $L$ (y por ende $L^*$) es invertible, resulta la definida positividad de $D$, lo que en este caso, por ser diagonal, equivale a $D>0$. ==Esto nos garantiza la existencia de la factorización de Cholesky==.
$$A=LD^{\frac{1}{2}}D^{\frac{1}{2}}L^*=CC^*$$
donde $C$ es triangular con elementos diagonales positivos.
### Base algorítmica de la factorización de Cholesky: Función $\text{chol}()^3$
La construcción de $C$ puede hacerse desde la factorización $LU$, sin embargo ese modo resulta más costoso y menos estable que la factorización de Cholesky. Antes, sin embargo, veremos como sería la base de un algoritmo (garantizando su ejecución sin interrupciones).
Asociando el algoritmo con una función $\text{chol}()^3$. Es decir $\text{chol}(A)$ devuelve el factor $C$.
Ya que $A=CC*$, con $C$ triangular inferior $c_{i,i}>0$, (de hecho, ya probamos la existencia de esta factorización). Intentaremos "despejar" $C$.
Escribimos
$$C=\begin{pmatrix}
	c_{1,1} & | & 0\\
	- & - & -\\
	C(2:n,1) & | & C(2:n,2:n)
\end{pmatrix}$$
donde $C(2:n,2:n)$ tiene las mismas características que $C$, sabiendo que
$$\begin{pmatrix}
	a_{1,1} & | & A(2:n,1)^*\\
	- & - & -\\
	A(2:n,1) & | & A(2:n,2:n)
\end{pmatrix} = $$
$$\begin{pmatrix}
	c_{1,1} & | & 0\\
	- & - & -\\
	C(2:n,1) & | & C(2:n,2:n)
\end{pmatrix}\
\begin{pmatrix}
	c_{1,1} & | & C(2:n,1)^*\\
	- & - & -\\
	0 & | & C(2:n,2:n)^*
\end{pmatrix}$$
se tiene, respectivamente, que los pasos a realizar son:
- **Paso 1**: $c_{1,1}^2=a_{1,1}\quad\rightarrow\quad c_{1,1}:=\sqrt{a_{1,1}}$
- **Paso 2**: $A(2:n,1)=c_{1,1}C(2:n,1)\quad\rightarrow\quad C(2:n,1):=A(2:n,1)/c_{1,1}$, lo cual nos permite "despejar" las incógnitas $C(1:n,1)$, puesto que $a_{1,1}>0$ por hipótesis, lo que permite elegir $c_{1,1}>0$. Ahora consideramos
- **Paso 3**: $A(2:n,2:n)=C(1:n,1)C(1:n,1)^*+C(2:n,2:n)C(2:n,2:n)^*\rightarrow$ $C(2:n,2:n):=\text{chol}(A(2:n,2:n)-C(1:n,1)C(1:n,1)^*)$ de donde el problema de resolver $\text{chol}(A(1:n,1:n))$ para una matriz de $n\times n$ se redujo a resolver $\text{chol}(A2:n,2:n)-C(2:n,1)C(2:n,1)^*)$ para una matriz de tamaño $(n-1)\times(n-1)$, lo que inductivamente (siempre que la matriz nueva se definida positiva) nos lleva al problema trivial de factorizar un número real positivo (matriz de $1\times1$ definida positiva) como el cuadrado de otro positivo.
>[!tip] Factorización de Cholesky en Python con `Numpy` y `Matlab`
>Con `Numpy` podemos obtener la factorización de Cholesky mediante `np.linal.cholesky()`, y en `Matlab` lo obtenemos mediante `chol()`.
### Propiedad de Reducción de Cholesky
>[!note] Propiedad de Reducción de Cholesky
>Si $A$ es $\text{SDP}$, entonces $A(2:n,2:n)-C(2:n,1)C(2:n,1)^*$, donde $C=(2:n,1)=A(2:n,1)/\sqrt{a_{1,1}}$, es $\text{DP}$.
#### Demostración
Claramente $A(2:n,2:n)-C(2:n,1)C(2:n,1)^*$ es Hermitiana (o simétrica, en el caso real). Por otro lado, y por hipótesis, para todo $0\neq u=(u_1,v)\in\mathbb{C}^n$ se tiene
$$0<u^*Au=\begin{pmatrix}u_1\\-\\v\end{pmatrix}^*\begin{pmatrix}a_{1,1} & | & A(2:n,1)^*\\-&-&-\\A(2:n,1)&|&A(2:n,2:n)\begin{pmatrix}u_1\\-\\v\end{pmatrix}\end{pmatrix}=$$
$$|u_1|^2a_{1,1}+u_1^*A(2:n,1)^*v+u_1v^*A(2:n,1)+v^*A(2:n,2:n)v$$
Elijamos $u_1=-\frac{A(2:n,1)^*v}{a_{1,1}}$ y resulta para todo $v\neq0$,
$$0<v^*(A(2:n,2:n)-C(2:n,1)C(2:n,1)^*)v$$
lo que completa la demostración.
### Implementación de factorización de Cholesky
En forma de algoritmo, Cholesky puede escribirse del siguiente modo
``` Plain Text
for j = i:n
	a_{j,j} := \sqrt{a_{j,j}}
		for i = j+1:n
			a_{i,j} := a_{i,j}/a_{j,j}
		for k = j+1:n
			for i = k:n
				a_{i,k} := a_{i,k}-a_{i,j}a_{k,j}
```
en donde se utiliza solo la información de la parte triangular inferior de $A$ (la parte superior es redundante) y en donde, además, se almacenan los datos de la matriz $C$.
# Ortogonalidad, Subespacios ortogonales y Proyecciones
## Ortogonalidad y subespacios ortogonales
### Vectores ortogonales
>[!note] Vectores ortogonales
>Dos vectores, $x$ e $y$, son ortogonales si el ángulo $\theta$ entre ellos es tal que $\cos{\theta}=0$ o, en términos del producto interno $\left<x,y\right>=0$. Esta relación la denotamos como $x\perp y$.
### Subespacio ortogonal
>[!note] Conjunto ortogonal
Sea $X=\{x_1,x_2,...\}$ un conjunto de vectores no nulos, este es un conjunto ortogonal si existen dos vectores ortogonales en el.
### Conjunto ortonormal
>[!note] Conjunto ortonormal
>Si el conjunto de vectores $X$ posee dos vectores ortogonales y cada vector es unitario, entonces es un conjunto ortonormal. Por lo tanto, para dos vectores cualesquiera en este conjunto ortonormal tenemos que $\left<x_i,x_j\right>=\delta_{ij}$, donde $\delta_{ij}=1$ si $i=j$ y $0$ en el resto de los casos.
>Además, decimos que un conjunto ortonormal es **completo** si no está contenido en ningún conjunto ortonormal mayor.
### Identidad Pitagórica
>[!note] Identidad Pitagórica
>Sean $x$ e $y$ dos vectores ortogonales en $\mathbb{R}^n$. Entonces,
>$$||x+y||^2=||x||^2+||y||^2$$

Podemos generalizarlo como, $\left|\left|\sum_{i=1}^kx_i\right|\right|^2=\sum_{i=1}^k||x_i||^2$ si $\{x_1,...,x_i\}$ es un conjunto ortogonal.
## Proyecciones ortogonales
### Expansiones ortogonales
Sea $X=\{x_1,...,x_n\}$ un conjunto ortogonal y $x$ un vector que pertenece al subespacio $X$. La ortogonalidad de los vectores en $X$ vuelve realmente fácil encontrar los coeficientes cuando $x$ es expresado como una combinación lineal de vectores en $X$. O sea, tenemos que
$$x=\alpha_1x_1+\alpha_2x_2+...+\alpha_nx_n=\sum_{i=1}^n\alpha_ix_i$$
Utilizando el producto interno con cualquier vector de $X$ tenemos que
$$\left<x,x_j\right>=\left<\sum_{i=1}^n\alpha_ix_i,x_j\right>=\sum_{i=1}^n\alpha_i\left<x_i,x_j\right>=\alpha_j\left<x_j,x_j\right>$$
puesto que $\left<x_i,x_j\right>=0$ para todo $i\neq j$. Esto muestra que $\alpha_j=\frac{\left<x,x_j\right>}{\left<x_j,x_j\right>}$ para cualquier $j=1,2,...,n$. Sustituyendo estos valores para los vectores de $X$ en la combinación lineal de estos tenemos que
$$x=\frac{\left<x,x_1\right>}{\left<x_1,x_1\right>}x_1+\frac{\left<x,x_2\right>}{\left<x_2,x_2\right>}x_2+\cdots+\frac{\left<x,x_n\right>}{\left<x_n,x_n\right>}x_n=\sum_{i=1}^n\frac{\left<x,x_i\right>}{\left<x_i,x_i\right>}x_i$$
Las cosas son aún más sencillas cuando $X=\{x_1,x_2...,x_n\}$ es un conjunto ortonormal. Entonces, tenemos que cada $\alpha_j=\left<x,x_j\right>$ y $x=\sum_{i=1}^n\left<x,x_i\right>x_i$.
### Función de proyección
Las expansiones ortogonales como la dada anteriormente pueden ser desarrolladas utilizando conceptos geométricos, en particular mediante la **función de proyección**.
>[!note] Función de proyección
>Sean $v$ y $x$ dos vectores en $\mathbb{R}^n$ tal que $v\neq0$. Entonces, la **proyección** de $x$ sobre el vector $v$ es dada por el vector función:
>$$\text{proj}_v(x)=\frac{\left<x,v\right>}{\left<v,v\right>}v,\quad\text{con }v\neq0$$

Geométricamente, se obtiene el **vector proyección** de un vector $x$ en la dirección de un vector $v$. Alternativamente, definimos a la función $\text{proj}_v(x)$ como la **proyección ortogonal** de $x$ en $v$. Cuando $v$ está normalizado, obtenemos que $\text{proj}_v(x)=\left<x,v\right>v=(v'x)v$.
>[!note] Componente
En análisis geométrico, la cantidad $\frac{v'x}{||v||}$ es usualmente referido como el **componente** de $x$ en la dirección de $v$ y es denotado como $\text{comp}_v(x)$. La longitud de la proyección de $x$ sobre $v$ es igual a el valor absoluto del componente de $x$ a lo largo de $v$:
$$||\text{proj}_v(x)||=\left(\frac{|v'x|}{v'v}\right)||v||=\frac{|v'x|}{||v||}=|\text{comp}_v(x)|$$
### Propiedades
Sean $v$ y $x$ dos vectores en $\mathbb{R}^n$. Asumimos $v\neq0$ y definimos:
$$P_v=\frac{vv'}{v'v}=\frac{1}{||v||^2}vv'$$
1. $P_v$ es una matriz de $\mathbb{R}^{n\times n}$ de rango uno que es simétrica e idempotente (igual a su cuadrado).
2. $\text{proj}_v(x)=P_vx$
3. $\text{proj}_v(x)$ es lineal en $x$
4. Si $x$ y $v$ son ortogonales, entonces $\text{proj}_v(x)=0$ y $P_vP_x=O$
5. $P_vv=v$
6. Si $x=\alpha v$, para algún $\alpha\in\mathbb{R}$, entonces $\text{proj}_v{x}=x$
7. Si $u=\alpha v$ para algún $\alpha\in\mathbb{R}-\{0\}$, entonces $P_u=P_v$ y $\text{proj}_v(x)=\text{proj}_u(x)$
### Proyección 
## Subespacios ortocomplementarios
## Teorema Fundamental del Álgebra Lineal
# Ortogonalidad en Métodos Directos
El método de eliminación se basa en triangulación por matrices triangulares (los multiplicadores). El crecimiento de los multiplicadores y de la matriz triangular superior resultante podría ser problemático. Por esa razón lo ideal sería realizar operaciones con matrices de norma controlada (idealmente de norma 1). Esos algoritmos son posibles y se basan en ideas de ortogonalidad.
## Factorización QR
$$A=aR,\qquad Ax=b$$
Dada una matriz $A\in\mathbb{C}^{m\times n}$, con columnas $\{a_1,a_2,...,a_n\}$ linealmente independientes ($m\geq n$), podemos considerar los subespacios anidados generados por sus columnas:
$$<a_1>\subset<a_1,a_2>\subset\cdots\subset<a_1,a_2,...,a_n>$$
En varias aplicaciones es de interés encontrar generadores **ortonormales** de esos subespacios. O sea, una colección $\{q_1,q_2,...,q_n\}$, tales que $q_i^*q_j=\delta_i^j$ y $<a_1,a_2,...,a_l>=<q_1,q_2,...,q_l>$, para todo $l$. Esto remite a **Gram-Schmidt**, el cual funcionaría de esta forma:
- $u_1=a_1\rightarrow q_1=\frac{u_1}{||u_1||_2}$
- $u_2=a_2-P_{u_1}(a_2)\rightarrow q_2 =\frac{u_2}{||u_2||_2}$
- $u_3=a_3-P_{u_1}(a_3)-P_{u_2}(a_3)\rightarrow q_3=\frac{u_3}{||u_3||_2}$
- $\cdots$
- $u_k=a_k-\sum_{i=1}^{k-1}P_{u_i}(a_k)\rightarrow q_k=\frac{u_k}{||u_k||_2}$
Matricialmente, y observando los subespacios anidados, es obvio que esto nos provee una matriz $R\in\mathbb{K}^{n\times n}$ triangular superior
$$R=\begin{pmatrix}
	r_{1,1} & r_{} & r_{} & \cdots & r_{} \\
	0 & r_{} & r_{} & \cdots & r_{}\\
	0 & 0 & r_{} &\cdots & r_{}\\
	\vdots & \vdots & \vdots & \vdots & \vdots \\
	0 & 0 & 0 & 0 & 0 
\end{pmatrix}$$
con $r_{i,i}>0,r_{i,i}\in\mathbb{R}$, tal que
$$\begin{pmatrix}
	    & | &     & | &        & | & \\
	    & | &     & | &        & | & \\
	a_1 & | & a_2 & | & \cdots & | & a_n\\
	    & | &     & | &        & | & \\
	    & | &     & | &        & | & \\ 
\end{pmatrix}=
\begin{pmatrix}
	    & | &     & | &        & | & \\
	    & | &     & | &        & | & \\
	q_1 & | & q_2 & | & \cdots & | & q_n\\
	    & | &     & | &        & | & \\
	    & | &     & | &        & | & \\
\end{pmatrix}
\begin{pmatrix}
	r_{1,1} & r_{} & r_{} & \cdots & r_{} \\
	0 & r_{} & r_{} & \cdots & r_{}\\
	0 & 0 & r_{} &\cdots & r_{}\\
	\vdots & \vdots & \vdots & \vdots & \vdots \\
	0 & 0 & 0 & 0 & 0 
\end{pmatrix}$$
de donde se observa que tanto los $q_i$ (elegidos ortonormales) como los $r_{i,j}$ pueden obtenerse desde las ecuaciones
$$\begin{cases}
	a_1 & = & r_{1,1}q_1\\
	a_2 & = & r_{1,2}q_1+r_{2,2}q_2\\
	\ &\ \vdots & \\
	a_n & = & r_{1,n}q_1+r_{2,n}q_2+\cdots+r_{n,n}q_n
\end{cases}$$
Notemos que en caso de que $m=n$, $Q$ resulta cuadrada ($Q\in\mathbb{K}^{n\times n}$), y unitaria (ortogonal en caso de ser una matriz real).
>[!note] Teorema Factorización QR Ortonormal-Triangular
>Dada $A\in\mathbb{K}^{m\times n},m\geq n$, con columnas l.i., es posible escribirla como producto de dos factores $Q\in\mathbb{K}^{m\times n}$ y $R\in\mathbb{K}^{n\times n}$. $Q$ con columnas ortonormales y $R$ triangular inferior.
### Observaciones
- Si utilizamos Gram-Schmidt vemos que $\mathbb{R}\ni r_{i,i}>0$.
- Consideremos la factorización $A=QR$: podemos multiplicar una columna arbitraria, digamos la columna $j$, de $Q$ por un complejo $z$ de la forma $e^{i\theta}\neq1$ y la fila $j$ por $z^{-1}=e^{-i\theta}$ obteniendo una nueva factorización $QR$. El resultado siguiente dice que salvo esa modificación, la factorización es única.
### Teorema de Factorización $QR$ única
>[!note] Teorema de Factorización $QR$ única con Diagonal Positiva
>Dada $A\in\mathbb{C}^{m\times n}$ con $m\geq n$ y rango máximo existe una única factorización $QR$ con la propiedad de que $\mathbb{R}\ni r_{i,i}>0$
### Algoritmo de Factorización QR
1. Uso **Gram-Schmidt** sobre cada una de las columnas de $A$, de esa forma obtengo vectores ortogonales entre sí, de norma 1.
2. Los vectores columnas obtenidos serán las columnas de $Q$, que tiene la misma dimensionalidad de $A$.
3. Como $Q$ es ortogonal (matriz con columnas ortonormales), $Q$ es inversible pero además $Q^t=Q^{-1}$.
4. Resuelvo la siguiente ecuación para hallar la matriz $R$
$$A=QR\implies Q^{-1}A=Q^{-1}QR\implies Q^tA=R$$