# Transformaciones lineales
Las **transformaciones lineales** son funciones entre $K$-espacios vectoriales que son compatibles con la estructura de estos espacios.
## Definiciones de Transformaciones lineales
### Definición formal
>[!note] Transformación lineal
>Sean $(V,+_V,\cdot_V)$ y $(W,+_W,\cdot_W)$ dos $K$-espacios vectoriales. Llamamos **transformación lineal** (u homomorfismo, o simplemente morfismo) a $f:V\to W$ si cumple:
>1. $f(v+_Vv')=f(v)+_Wf(v'),\forall v,v'\in V$
>2. $f(\lambda\cdot_V v)=\lambda_Wf(v),\forall\lambda\in K,\forall v\in V$

O sea, tenemos que <mark style="background: #FFF3A3A6;">una aplicación es una transformación lineal si "preserva" las dos operaciones básicas de un espacio vectorial, la suma vectorial y el producto por un escalar.</mark> De este modo, obteniendo una mayor generalidad, para escalares cualesquiera $a_i\in K$ y vectores cualesquiera $v_i\in V$ llegamos a la <mark style="background: #FFB86CA6;">propiedad básica de las aplicaciones lineales</mark>:
$$F(a_1v_1+a_2v_2+\cdots+a_nv_n)=a_1F(v_1)+a_2F(v_2)+\cdots+a_nF(v_n)$$
#### Observación: <mark style="background: #FF5582A6;">"3er condición" de una TL</mark>
Si $f:V\to W$ es una transformación lineal, entonces $f(0_V)=0_W$.
##### Demostración
Puesto que $f(0_V)=f(0_V+0_V)=f(0_V)+f(0_V)$, entonces
$$0_W=f(0_V)+(-f(0_V))=\left(f(0_V)+f(0_V)\right)+(-f(0_V))$$
$$=f(0_V)+(f(0_V)+(-f(0_V)))=f(0_V)+0_W=f(0_V)$$
#### Ejemplo
1. Sean $V$ y $W$ dos $K$-espacios vectoriales. Entonces $0:V\to W$, definida por $0(x)=0_W,\forall x\in V$, es una transformación lineal.
2. Si $V$ es un $K$-espacio vectorial, $id:V\to V$ definida por $id(x)=x$ es una transformación lineal.
3. Sea $A\in K^{m\times n}$. Entonces $f_A:K^n\to K^m$ definida por $f_A(x)=(A\cdot x^t)^t$ es una transformación lineal.
4. $f:K[X]\to K[X],f(P)=P'$ es una transformación lineal.
5. $F:C(\mathbb{R})\to\mathbb{R}$, donde $C(\mathbb{R})=\{f:\mathbb{R}\to\mathbb{R}:f$ es continua$\}, F(g)=\int_0^1g(x)\ dx$ es una transformación lineal.
### Definición matricial
>[!note] Definición matricial de una transformación lineal
>Una matriz $A\in\mathbb{K}^{m\times n}$ define una función $T_A:\mathbb{R}^n\to\mathbb{R}^m$ dada por
>$$T_A(v)=Av$$
>siendo $v$ una matriz columna que representa a un vector.
#### Ejemplo
Si $A=\begin{pmatrix} 1 & 3 & 0 \\ 2 & -1 & -4 \end{pmatrix}\in\mathbb{R}^{2\times3}$ y $v=(1,0,-2)$, entonces
$$T_A(v)=\begin{pmatrix} 1 & 3 & 0 \\ 2 & -1 & -4 \end{pmatrix}\begin{pmatrix} 1 \\ 0 \\ -2 \end{pmatrix}=\begin{pmatrix} 1 \\ 10 \end{pmatrix}$$
Si multiplicamos $A$ por un vector genérico $v(x_1,x_2,x_3)\in\mathbb{R}^3$, obtenemos
$$T_A(v)=\begin{pmatrix} 1 & 3 & 0 \\ 2 & -1 & -4 \end{pmatrix}\begin{pmatrix} x_1 \\ x_2 \\ x_3\end{pmatrix}=\begin{pmatrix} x_1+3x_2+0x_3 \\ 2x_1-x_2-4x_3 \end{pmatrix}$$
Entonces podemos escribir
$$T_A(x_1,x_2,x_3)=(x_1+3_x2,2x_1-x_2-4x_3)$$
#### Ejercicio
Escribir la función $f(x_1,x_2)=(3x_1,x_2-x_1,x_2+2x_1)$ de $\mathbb{R}^2$ en $\mathbb{R}^3$ en forma matricial.
$$T_A(x_1,x_2)=\begin{pmatrix} 3x_1 \\ x_2-x_1 \\ x_2+2x_1 \end{pmatrix}$$
Puesto que esta expresión matricial se obtiene a partir del producto $Av,v=(x_1,x_2)^t$, tenemos que la matriz $A$ que define a la función $f$ es la siguiente:
$$A=\begin{pmatrix}
	3 & 0\\
	-1 & 1\\
	2 & 1
\end{pmatrix}$$
## Propiedades
### Preservación de subespacios

^68c3d8

>[!note] Preservación de subespacios
>Sea $f:V\to W$ una transformación lineal. Entonces:
>1. **Preservación de subespacios bajo la imagen**: Si $S$ es un subespacio de $V$, entonces $f(S)$ es un subespacio de $W$.
>2. **Preservación de subespacios bajo la pre-imagen**: Si $T$ es un subespacio de $W$, entonces $f^{-1}(W)$ es un subespacio de $V$.
#### Demostración
1. Sea $S\subset V$ un subespacio y consideremos $f(S)=\{w\in W:\exists s\in S, f(s)=w\}$.
	1. $0_W\in f(S)$, puesto que $f(0_V)=0_W$ y $0_V\in S$.
	2. Sean $w,w'\in f(S)$. Entonces existen $s,s'\in S$ tales que $w=f(s)$ y $w'=f(s')$. Luego $w+w'=f(s)+f(s')=f(s+s')\in f(S)$, puesto que $s+s'\in S$.
	3. Sean $\lambda\in K$ y $w\in f(S)$. Existe $s\in S$ tal que $w=f(s)$. Entonces $\lambda\cdot w=\lambda f(s)=f(\lambda\cdot s)\in f(S)$, puesto que $\lambda\cdot s\in S$.
2. Sea $T$ un subespacio de $W$ y consideremos $f^{-1}(T)=\{v\in V: f(v)\in T\}$
	1. $0_V\in f^{-1}(T)$, puesto que $f(0_V)=0_W\in T$.
	2. Sean $v,v'\in f^{-1}(T)$. Entonces $f(v),f(v')\in T$ y, por lo tanto, $f(v+v')=f(v)+f(v')\in T$. Luego $v+v'\in f^{-1}(T)$.
	3. Sean $\lambda\in K, v\in f^{-1}(T)$. Entonces $f(v)\in T$ y, en consecuencia, $f(\lambda\cdot v)=\lambda\cdot f(v)\in T$. Luego $\lambda\cdot v\in f^{-1}(T)$
### Consecuencia de la definición formal
<mark style="background: #ADCCFFA6;">Una transformación lineal preserva combinaciones lineales.</mark> Debido a esto, una transformación queda unívocamente determinada por los valores que toma en los elementos de una base cualquiera de su dominio.
#### Ejemplo
Hallar, si es posible, una transformación lineal $f:\mathbb{R}^2\to\mathbb{R}^2$ que verifique $f(1,1)=(0,1)$ y $f(1,0)=(2,3)$.
Dado $(x_1,x_2)\in\mathbb{R}^2$ se tiene que $(x_1,x_2)=x_2(1,1)+(x_1-x_2)(1,0)$. Entonces, si $f$ verifica lo pedido, debe ser
$$f(x_1,x_2)=x_2\cdot f(1,1)+(x_1-x_2)\cdot f(1,0)=x_2\cdot (0,1)+(x_1-x_2)\cdot(2,3)$$
$$=(2x_1-2x_2,3x_1-2x_2)$$
Además, es fácil ver que esta función es una transformación lineal que vale $f(1,1)=(0,1)$ y $f(1,0)=(2,3)$.
Luego, $f(x_1,x_2)=(2x_1-2x_2,3x_1-2x_2)$ es la única transformación lineal que satisface lo pedido.
### Unicidad de la transformación lineal
>[!note] Unicidad de la transformación lineal
>Sean $V$ y $W$ dos $K$-espacios vectoriales, $V$ de dimensión finita. Sea $B=\{v_1,...,v_n\}$ una base de $V$ y sean $w_1,...,w_n\in W$ *vectores arbitrarios*. Entonces **existe una única transformación lineal** $f:V\to W$ tal que $f(v_i)=w_i$ para cada $1\leq i\leq n$.

El hecho de que los vectores $w_1,...,w_n$ son completamente arbitrarios radica en el hecho de que pueden ser linealmente dependientes o incluso iguales entre sí.
O sea, <mark style="background: #FFB8EBA6;">si una transformación lineal dada es definida a partir de su aplicación sobre un conjunto de vectores que constituyen una base del dominio, sin importar el resultado de cada aplicación (denominados vectores arbitrarios en definición), tenemos que dicha TL es única.</mark>
#### Demostración
##### Existencia
Dado $v\in V$ existen únicos $\alpha_1,...,\alpha_n\in K$ tales que $v=\sum_{i=1}^n\alpha_iv_i$, es decir, $(\alpha_1,...,\alpha_n)=(v)_B$ es el vector de coordenadas de $v$ en la base $B$. Definimos
$$f(v)=\sum_{i=1}^n\alpha_iw_i$$
Verificamos que $f$ es una transformación lineal. Sean $v,v'\in V$. Supongamos que $v=\sum_{i=1}^n\alpha_iv_i$ y $v'=\sum_{i=1}^n\alpha'_iv_i$. Entonces
$$v+v'=\sum_{i=1}^n\alpha_iv_i+\sum_{i=1}^n\alpha'_iv_i=\sum_{i=1}^n(\alpha_i+\alpha'_i)v_i$$
y, en consecuencia,
$$f(v+v')=\sum_{i=1}^n(\alpha_i+\alpha'_i)w_i=\sum_{i=1}^n\alpha_iw_i+\sum_{i=1}^n\alpha'_iw_i=f(v)+f(v')$$
De manera análoga se prueba que $f(\lambda v)=\lambda f(v), \forall \lambda\in K,\forall v\in V$.
##### Unicidad
Supongamos que $f$ y $g$ son dos transformaciones lineales de $V$ en $W$ tales que $f(v_i)=w_i$ y $g(v_i)=w_i$ para cada $1\leq i\leq n$. Entonces, dado $v\in V$, si $v=\sum_{i=1}^n\alpha_iv_i$, por la linealidad de $f$ y $g$ se tiene que
$$f(v)=\sum_{i=1}^n\alpha_if(v_i)=\sum_{i=1}^n\alpha_ig(v_i)=g(v)$$
Luego, $f(v)=g(v)$ para todo $v\in V$, de donde $f=g$.
##### Observación
Con una demostración análoga, de verifica que, si $V$ y $W$ son dos $K$-espacios vectoriales ($V$ no necesariamente de dimensión finita), $B=\{v_i:i\in I\}$ una base de $V$ y $\{w_i:i\in I\}\subset W$, existe una única transformación lineal $f:V\to W$ tal que $f(v_i)=w_i,\forall i\in I$.
## Clasificación de transformaciones lineales
>[!note] Clasificación funcional de transformaciones lineales
>Sean $V$ y $W$ dos $K$-espacios vectoriales, y sea $f:V\to W$ una transformación lineal. Se dice que:
>- $f$ es un *monomorfismo* si $f$ es **inyectiva**
>- $f$ es un *epimorfismo* si $f$ es **suryectiva**
>- $f$ es *isomorfismo* si $f$ es **biyectiva**
### Endomorfismo y automorfismo
>[!note] Endomorfismo y automorfismo
>Sea $V$ un $K$-espacio vectorial. Una transformación lineal $f:V\to V$ se llama **endomorfismo** de $V$. Si $f$ es un endomorfismo que es además un isomorfismo, entonces se dice que es un **automorfismo**.
## Núcleo e imagen de una Transformación Lineal
### Núcleo
>[!note] Núcleo
>Sean $V$ y $W$ dos $K$-espacios vectoriales, y sea $f:V\to W$ una transformación lineal. Se llama **núcleo** de $f$ al conjunto $\text{Nu}(f)=\{v\in V: f(v)=0\}=f^{-1}(\{0\})$
#### Importancia
A una transformación lineal $f:V\to W$ podemos asociarle un subespacio de $V$, llamado su <mark style="background: #ADCCFFA6;">núcleo, que de alguna manera mide el tamaño de la pre-imagen</mark> por $f$ de un elemento de su imagen. En particular, conocer este subespacio nos permitirá determinar si $f$ es inyectiva.
#### Nulidad
>[!note] Nulidad
>Sea $F:V\to U$ una transformación lineal, se define a la nulidad de $F$ como la <mark style="background: #FFF3A3A6;">dimensión de su núcleo</mark>. O sea,
>$$\text{nulidad }F=\text{dim}(\text{Nu }F)$$
#### Observación
Si $f:V\to W$ es una transformación lineal, $\text{Nu}(f)$ es un subespacio de $V$, puesto que es la pre-imagen por $f$ del subespacio $\{0\}\subset W$.
#### Ejemplo
Sea $f:\mathbb{R}^3\to\mathbb{R}^2$ la TL definida por $f(x_1,x_2,x_3)=(x_1,x_2)$. Entonces
$$\text{Nu}(f)=\{(x_1,x_2,x_3)\in\mathbb{R}^3:f(x_1,x_2,x_3)=0\}$$
$$=\{(x_1,x_2,x_3)\in\mathbb{R}^3:x_1=x_2=0\}$$
$$=\left<(0,0,1)\right>$$
#### Núcleo de un monomorfismo
>[!note] Condición del núcleo de un monomorfismo
>Sea $f:V\to W$ una transformación lineal. Entonces
>$$f\ \text{es monomorfismo}\iff\text{Nu}(f)=\{0\}$$
##### Demostración
###### $\implies$
Si $f$ es un monomorfismo, entonces es una función inyectiva. En particular, <mark style="background: #D2B3FFA6;">existe a lo sumo un elemento del dominio tal que la transformación nos da al vector nulo</mark>. Puesto que $f(0)=0$, debe ser $v=0$. Luego, $\text{Nu}(f)=\{0\}$.
###### $\Longleftarrow$
Sean $v,v'\in V$ Supongamos que $f(v)=f(v')$. Entonces $f(v-v')=f(v)-f(v')=0$, con lo que $v-v'\in\text{Nu}(f)$ y por lo tanto, la hipótesis $\text{Nu}(f)=\{0\}$ implica que $v-v'=0$, es decir, $v=v'$. Luego $f$ es inyectiva.
### Imagen
>[!note] Imagen
>Sea $f:V\to W$ una transformación lineal. Su **imagen** se define como
>$$\text{Im}(f)=\{w\in W:\exists v\in V,f(v)=w\}$$

Por la [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#^68c3d8|propiedad de preservación de subespacios]], tenemos que la imagen de una TL $f:V\to W$ resulta ser un subespacio de $W$.
#### Rango
>[!note] Rango
>Sea $F:V\to U$ una transformación lineal. Se define el rango de $F$ como la <mark style="background: #FFF3A3A6;">dimensión de su imagen</mark>. O sea
>$$\text{rango }F=\text{dim}(\text{Im }F)$$
#### Ejemplo 1: cálculo de imagen por definición
Hallar la imagen de la transformación lineal $f:\mathbb{R}^3\to\mathbb{R}^3$ definida como $f(x_1,x_2,x_3)=(x_1-x_2,-x_1+x_2,2x_1-2x_2+x_3)$.
Por definición,
$$\text{Im}(f)=\{y\in\mathbb{R}^3:\exists x\in\mathbb{R}^3,f(x)=y\}$$
$$=\{y\in\mathbb{R}^3:\exists(x_1,x_2,x_3)\in\mathbb{R}^3,(x_1-x_2,-x_1+x_2,2x_1-2x_2+x_3)=y\}$$
Entonces, un elemento de $y$ pertenece a $\text{Im}(f)$ si y sólo si es de la forma:
$$y=(x_1-x_2,-x_1+x_2,2x_1-2x_2+x_3)$$
$$=(x_1,-x_1,2x_1)+(-x_2,x_2,-2x_2)+(0,0,x_3)$$
$$=x_1\cdot(1,-1,2)+x_2\cdot(-1,1,-2)+x_3\cdot(0,0,1)$$
Luego
$$\text{Im}(f)=<(1,-1,2),(-1,1,-2),(0,0,1)>=<(1,-1,2),(0,0,1)>$$
#### Ejemplo 2: cálculo de imagen por canónicos
Otra manera de calcular la imagen de $f$, teniendo en cuenta que es una transformación lineal, es la siguiente:
Considerando un sistema de generadores de $\mathbb{R}^3$, por ejemplo la base canónica $\{e_1,e_2,e_3\}$. Para cada $x\in\mathbb{R}^3$ se tiene que $x=x_1\cdot e_1+x_2\cdot e_2+x_3\cdot e_3$, de donde resulta que
$$f(x)=x_1\cdot f(e_1)+x_2\cdot f(e_2)+x_3\cdot f(e_3)$$
Luego,
$$\text{Im}(f)=\{f(x):x\in\mathbb{R}^3\}=<f(e_1),f(e_2),f(e_3)>$$
$$=<(1,-1,2),(-1,1,-2),(0,0,1)>=<(1,-1,2),(0,0,1)>$$
#### Propiedades
##### Transformación de los Generadores
>[!note] Propiedad de Transformación de los Generadores de $V$
>Sean $V$ y $W$ dos $K$-espacios vectoriales y sea $f:V\to W$ una transformación lineal. Entonces, si $\{v_i:i\in I\}$ es un sistema de generadores de $V,\{f(v_i):i\in I\}$ es un sistema de generadores de $\text{Im}(f)$.
###### Demostración
Por definición, $\text{Im}(f)=\{w\in W:\exists v\in V, f(v)=w\}=\{f(v):v\in V\}$.
Si $\{v_i:i\in I\}$ es un sistema de generadores de $V$, para cada $v\in V$, existen $i_1,...,i_n\in I$ y elementos $\alpha_{i_j}\in K$ tales que $v=\sum_{j=1}^n\alpha_{i_j}v_{i_j}$. Luego
$$f(v)=\sum_{j=1}^n\alpha_{i_j}f(v_{i_j})\in\left<\{f(v_i):i\in I\}\right>$$
Esto prueba que $\text{Im}(f)\subset\left<\{f(v_i):i\in I\}\right>$. Es claro que vale la otra inclusión, ya que $f(v_i)\in\text{Im}(f)$ para cada $i\in I$.
Luego, $\text{Im}(f)=\left<\{f(v_i):i\in I\}\right>$.
###### Corolario de Rango y Nulidad
>[!info] Corolario de Rango y Nulidad
>Sean $V$ y $W$ dos $K$-espacios vectoriales y sea $f:V\to W$ una TL. Si $V$ es de dimensión finita, entonces $\text{Im}(f)$ también lo es y se tiene que $\text{dim}(\text{Im}(f))\leq\text{dim}V$.
###### Consecuencia para un epimorfismo
>[!info] Consecuencia de la transformación de los Generadores para un epimorfismo
>Si $f:V\to W$ es un **epimorfismo**, y $\{v_i:i\in I\}$ es un sistema de generadores de $V$, entonces $\{f(v_i):i\in I\}$ es un sistema de generadores de $W$.

^ccd9dc

###### Ejemplo
Sean $S$ y $T$ los subespacios de $\mathbb{R}^3$ definidos por $S=\{x\in\mathbb{R}^3:x_1-x_2=0\}$ y $T=\{x_in\mathbb{R}^3:x_3=0\}$. Hallar una TL $f:\mathbb{R}^3\to\mathbb{R}^3$ tal que $f(S)=T$.
Para definir una TL $f:\mathbb{R}^3\to\mathbb{R}^3$ basta con especificar los valores que toma sobre los elementos de una base de $\mathbb{R}^3$.
Por lo tanto, consideramos la base de $S$, por ejemplo $\{(1,1,0),(0,0,1)\}$, y la extendemos a una base de $\mathbb{R}^3$, por ejemplo $\{(1,1,0),(0,0,1),(1,0,0)\}$. Teniendo en cuenta que $T=<(1,0,0)(0,1,0)>$, definimos:
$$f(1,1,0)=(1,0,0),\qquad f(0,0,1)=(0,1,0),\qquad f(1,0,0)=(0,0,1)$$
Entonces
$$f(S)=<f(1,1,0)=(1,0,0),\quad f(0,0,1)=(0,1,0),\quad f(1,0,0)=(0,0,1)$$

>[!warning] La transformación de los generadores no necesariamente es una base de $Im(f)$
>Observemos que si $f:V\to W$ es un epimorfismo y $\{v_i:i\in I\}$ es una base de $V$, entonces $\{f(v_i):i\in I\}$ no es necesariamente una base de $\text{Im}(f)$. Por el corolario anterior, es un sistema de generadores, pero podría no ser un conjunto linealmente independiente. Una transformación lineal arbitraria no preserva independencia lineal.
##### Preservación de la independencia lineal bajo monomorfismos
>[!note] Preservación de la Independencia Lineal bajo monomorfismos
>Sean $V$ y $W$ dos $K$-espacios vectoriales y sea $f:V\to W$ un **monomorfismo**. Entonces, si $\{v_i:i\in I\}\subset V$ es un conjunto linealmente independiente, $\{f(v_i):i\in I\}\subset W$ es un conjunto linealmente independiente.
###### Demostración
Supongamos que una combinación lineal de $\{f(v_i):i\in I\}$ satisface $\sum_{i\in I}\alpha_if(v_i)=0$. Como $f$ es una TL, entonces $f\left(\sum_{i\in I}\alpha_iv_i\right)=0$, y como es un monomorfismo, debe ser $\sum_{i\in I}\alpha_iv_i=0$. La independencia lineal de $\{v_i:i\in I\}$ implica que $\alpha_i=0,\forall i\in I$.
###### Corolario de preservación de la base bajo monomorfismos
>[!info] Preservación de la Base bajo monomorfismos
>Si $f:V\to W$ es un **monomorfismo** y $B=\{v_i:i\in I\}$ es una base de $V$, entonces $\{f(v_i):i\in I\}$ es una base de $\text{Im}(f)$. En particular, si $V$ es un $K$-espacio vectorial de dimensión finita, $\text{dim}(\text{Im}(f))=\text{dim}V$.

^9f5e2c

###### Corolario de preservación de la Base y Dimensión bajo Isomorfismos
Considerando que un isomorfismo es una TL que es a la vez epimorfismo y monomorfismo, conociendo la [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#^ccd9dc|consecuencia de la transformación de los Generadores para un epimorfismo]] y la [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#^9f5e2c|preservación de la base bajo monomorfismos]], podemos deducir que:
>[!info] Preservación de la Base y Dimensión bajo Isomorfismos
>Sean $V$ y $W$ dos $K$-espacios vectoriales y sea $f:V\to W$ un isomorfismo. Entonces para toda base $B$ de $V$, $f(B)$ es una base de $W$. En particular, si $V$ es de dimensión finita, $W$ también lo es y $\text{dim}V=\text{dim}W$.
### Espacios vectoriales de dimensión finita
Al estudiar espacios vectoriales de dimensión finita, podemos trabajar en un $K$-espacio vectorial arbitrario de dimensión $n$ "como si fuese" $K^n$ simplemente considerando vectores de coordenadas. La noción de isomorfismo nos permite formalizar esta idea.
#### Isomorfismo canónico de un espacio vectorial
>[!note] Isomorfismo canónico de un espacio vectorial
>Sea $V$ un $K$-espacio vectorial de dimensión $n$. Entonces existe un isomorfismo $f:V\to K^n$
##### Demostración
Sea $B=\{v_1,...,v_n\}$ una base de $V$. Dado $x\in V$, existe únicos $x_1,...,x_n\in K$ tales que $x=\sum_{i=1}^nx_iv_i$. Definimos:
$$f:V\to K^n,f(x)=(x_1,...,x_n)$$
Sean $x,y\in V$. Si $x=\sum_{i=1}^nx_iv_i$ e $y=\sum_{i=1}^ny_iv_i$, entonces $x+y=\sum_{i=1}^n(x_i+y_i)v_i$. Luego
$$f(x+y)=(x_1+y_1,...,x_n+y_n)=(x_1,...,x_n)+(y_1,...,y_n)=f(x)+f(y)$$
En forma análoga, se prueba que $f(\lambda\cdot x)=\lambda\cdot f(x)$ para cada $\lambda\in K$ y cada $x\in V$.
- Es claro que si $f(x)=0$, entonces $x=0$, de donde $f$ es un monomorfismo.
- Finalmente, dado $(x_1,...,x_n)\in K^n$, consideramos $x=\sum_{i=1}^nx_iv_i\in V$. Se tiene que $f(x)=(x_1,...,x_n)$. Luego, $f$ es un epimorfismo.
En consecuencia, $f$ es un isomorfismo.
#### Ejemplo
Sea $K_3[X]=\{P\in K[X]:P=0\lor\text{gr}(P)\leq3\}$, que es un $K$-espacio vectorial de dimensión 4.
Un isomorfismo $f:K_3[X]\to K^4$ puede definirse como sigue:
Si $P=\sum_{i=0}^3a_iX^i$, entonces $f(P)=(a_0,a_1,a_2,a_3)$, lo que corresponde a considerar en la demostración anterior la base $B=\{1,X,X^2,X^3\}$ de $K_3[X]$.
#### Observación
Observar que, considerando la aplicación $f$ definida en la demostración anterior, es tomar coordenadas en la base $B$, lo que nos permite trabajar con coordenadas de una base en el siguiente sentido:
1. $\{w_1,...,w_s\}$ es l.i. en $V\iff\{f(w_1),...,f(w_s)\}$ es l.i. de $K^n$.
2. $\{w_1,...,w_r\}$ es un sistema de generadores de $V\iff\{f(w_1),...,f(w_r)\}$ es un sistema de generadores de $K^n$.
3. $\{w_1,...,w_n\}$ es una base de $V\iff\{f(w_1),...,f(w_n)\}$ es una base de $K^n$.
##### Ejemplo
Para decidir si $\{X^2-X+1,X^2-3\cdot X+5,2\cdot X^2+2\cdot X-3\}$ es una base de $\mathbb{R}_2[X]$, bastará ver si $\{(1,-1,1),(1,-3,5),(2,2,-3)\}$ es una base de $\mathbb{R}^3$ para lo que se puede usar el método de triangulación.
### Teorema de Dimensión
Este teorema relaciona las dimensiones del núcleo y de la imagen de una TL con la de su dominio.
>[!note] Teorema de la dimensión para transformaciones lineales
>Sean $V$ y $W$ dos $K$-espacios vectoriales, $V$ de dimensión finita, y sea $f:V\to W$ una transformación lineal. Entonces
>$$\text{dim }V=\text{rango }F+\text{nulidad }F$$
#### Demostración
Sean $n=\text{dim}V$ y $r=\text{dim}(\text{Nu}(f))$.
Si $r=n$, entonces $f\equiv0$ y $\text{dim}(\text{Im}(f))=0$. Por lo tanto, el teorema vale.
Si $r=0$, entonces $f$ es un monomorfismo. En este caso, si $B$ es una base de $V$, entonces $f(B)$ es una base de $\text{Im}(f)$. Luego $\text{dim}(\text{Im}(f))=\text{dim}V$ ([[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#^9f5e2c|preservación de la base bajo monomorfismos]]), y el teorema vale.
Supongamos ahora que $0<r<n$. Sea $\{v_1,...,v_r\}$ una base de $\text{Nu}(f)$. Sean $v_{r+1},...,v_n$ en $V$ tales que $\{v_1,...,v_r,v_{r+1},...,v_n\}$ es una base de $V$. Verificamos entonces que $\{f(v_{r+1},...,f{v_n})\}$ es una base de $\text{Im}(f)$, de donde se deduce inmediatamente el teorema:
- Puesto que $\{v_1,...,v_r,v_{r+1},...,v_n\}$ es una base de $V$, se tiene que$$\text{Im}(f)=<f(v_1),...,f(v_r),f(v_{r+1}),...,f(v_n)>=f(v_{r+1},...,f(v_n)>,$$pues $f(v_i)=0$ para $1\leq i\leq r$.
- Sean $\alpha_{r+1},...,\alpha_n\in K$ tales que $\sum_{i=r+1}^n\alpha_if(v_i)=0$. Entonces $f\left(\sum_{i=r+1}^n\alpha_iv_i\right)=0$, es decir, $\sum_{i=r+1}^n\alpha_iv_i\in\text{Nu}(f)$. Como $\{v_1,...,v_r\}$ es una base de $\text{Nu}(f)$, existen $\alpha_1,...,\alpha_r\in K$ tales que $$\sum_{i=r+1}^n\alpha_iv_i=\sum_{i=r+1}^n\alpha_iv_i\iff\sum_{i=r+1}^n(-\alpha_i)v_i+\sum_{i=r+1}^n\alpha_iv_i=0$$Como $\{v_1,...,v_n\}$ es un conjunto linealmente independiente, $\alpha_i=0$ para cada $1\leq i\leq n$. En particular, $\alpha_i=0$ para $i=r+1,...,n$. Luego, $\{f(v_{r+1}),...,f(v_n)\}$ es un conjunto linealmente independiente.
De este modo se demuestra que si una TL entre dos espacios vectoriales de dimensión $n$ es inyectiva (respectivamente suryectiva), entonces también es suryectiva (respectivamente inyectiva).
#### Corolario
>[!info] Equivalencia de morfismos para transformaciones lineales
>Sean $V$ y $W$ dos $K$-espacios vectoriales de dimensión $n$ y sea $f:V\to W$ una transformación lineal. Son equivalentes:
>1. $f$ es un isomorfismo
>2. $f$ es un monomorfismo
>3. $f$ es un epimorfismo
##### Demostración
###### $1.\implies2.$
Por definición
###### $2.\implies3.$
Por el [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#Teorema de Dimensión|teorema de dimensión]], $n=\text{dim}V=\text{dim}(\text{Nu}(f))+\text{dim}(\text{Im}(f))$, y como $f$ es un monomorfismo, $\text{dim}(\text{Nu}(f))=0$. Entonces $\text{dim}(\text{Im}(f))=n=\text{dim}W$, de donde $\text{Im}(f)=W$.
###### $3.\implies1.$
Por el [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#Teorema de Dimensión|teorema de dimensión]], y considerando que $f$ es epimorfismo, se tiene que $n=\text{dim}V=\text{dim}(\text{Nu}(f))+\text{dim}(\text{Im}(f))=\text{dim}(\text{Nu}(f))+n$. Esto implica que $\text{dim}(\text{Nu}(f))=0$, con lo cual, $\text{Nu}(f)=\{0\}$ y $f$ es un monomorfismo. Siendo epimorfismo y monomorfismo, resulta que $f$ es un isomorfismo.
###### Advertencia
>[!warning] Corolario válido sólo para espacios de dimensión finita
>Como estos resultados son consecuencias del [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#Teorema de Dimensión|teorema de dimensión]], es esencial la hipótesis de que los espacios vectoriales sean de dimensión finita. Estos resultados no valen para TL definidas en espacios de dimensión infinita.
#### Ejemplo
Sea $V=K[X]$
1. Sea $f:K[X]\to K[X],f(P)=P'$, que es una TL.
	- $f$ es epimorfismo: Sea $Q=\sum_{i=0}^na_iX^i$. Entonces $f\left(\sum_{i=1}^{n+1}\frac{a_i}{i}X^i\right)=Q$.
	- $f$ no es monomorfismo: $f(1)=0$, pero $1\neq0$.
2. Sea $g:K[X]\to K[X],g(P)=X\cdot P$.
	- $g$ es monomorfismo: Si $f(P)=X\cdot P=0$, entonces $P=0$.
	- $g$ no es epimorfismo: $1\not\in\text{Im}(f)$. 
### Co-núcleo y Co-imagen
#### Co-núcleo
>[!note] Co-núcleo
>Espacio cociente del codominio con respecto a la imagen de la transformación. De este modo, sea una transformación lineal $T:V\to W$, el co-núcleo es el espacio cociente:$$\text{conúcleo}(T)=W\setminus\text{Im}(T)$$

O sea, el co-núcleo mide cuántos elementos del codominio no están cubiertos por la imagen de $T$.
#### Co-imagen
>[!note] Co-imagen
>Espacio de partida de la transformación, mapeado a través de la transformación $T$ (la imagen de la transformación). Sea una transformación lineal $T:V\to W$, la co-imagen es:$$\text{coimagen}(T)=\text{Im}(T)$$

O sea, la co-imagen representa todos los resultados posibles al aplicar $T$ a los elementos del dominio (i.e., la imagen de $T$).
## Composición de Transformaciones Lineales
La composición de funciones usual puede realizarse, en particular, entre dos TL, generando una nueva transformación lineal.
### Definición
>[!note] Composición de transformaciones lineales
>Sean $V$, $W$ y $Z$ $K$-espacios vectoriales. Sean $f:V\to W$ y $g:W\to Z$ TLs. Entonces $g\circ f:V\to Z$ es una TL.

Como la composición de TLs genera una TL, esta respeta la suma de vectores y la multiplicación por escalares. Entonces, sean $v,v'\in V$ tenemos que
$$g\circ f(v+v')=g(f(v+v'))=g(f(v)+f(v'))$$
$$=g(f(v))+g(f(v'))=g\circ f(v)+g\circ f(v')$$
Y también, si $\lambda\in K$ y $v\in V$, se tiene que
$$g\circ f(\lambda\cdot v)=g(f(\lambda\cdot v))=g(\lambda\cdot f(v))=\lambda\cdot g(f(v))=\lambda\cdot(g\circ f(v))$$
### Matriz de la composición de transformaciones lineales
La composición de dos TLs "se traduce" como la multiplicación de sus matrices.
>[!note] Matriz de una composición de transformaciones lineales
>Sean $V$, $W$ y $U$ tres $K$-espacios vectoriales de dimensión finita. Sean $B_1$, $B_2$ y $B_3$ bases de $V$, $W$ y $U$ respectivamente. Sean $f:V\to W$ y $g:W\to U$ transformaciones lineales. Entonces
>$$|g\circ f|_{B_1B_3}=|g|_{B_2B_3}\cdot|f|_{B_1B_2}$$
#### Demostración
Sean $n=\text{dim}V$, $m=\text{dim}W$ y $r=\text{dim}U$. Entonces $|g|_{B_2B_3}\in K^{r\times m}$ y $|f|_{B_1B_2}\in K^{m\times n}$, con lo que $|g|_{B_2B_3}\cdot|f|_{B_1B_2}\in K^{m\times n}$. Además, $|g\circ f|_{B_1B_3}\in K^{r\times n}$. Para cada $x\in V$ se tiene que
$$|g|_{B_2B_3}\cdot|f|_{B_1B_2}\cdot(x)_{B_1}=|g|_{B_2B_3}(f(x))_{B_2}=g(f(x))_{B_3}$$
$$=(g\circ f(x))_{B_3}=|g\circ f|_{B_1B_3}\cdot(x)_{B_1}$$
Luego, $|g|_{B_2B_3}\cdot|f|_{B_1B_2}=|g\circ f|_{B_1B_3}$.
#### Inversión de isomorfismos en bases de espacios vectoriales
>[!info] Inversión de isomorfismos en Bases de Espacios vectoriales
>Sean $V$ y $W$ dos $K$-espacios vectoriales de dimensión finita, y sean $B_1$ y $B_2$ bases de $V$ y $W$ respectivamente. Sea $f:V\to W$ un isomorfismo. Entonces $|f^{-1}|_{B_2B_1}=(|f|_{B_1B_2})^{-1}$.
##### Demostración
Se deduce inmediatamente aplicando la proposición anterior a $f^{-1}\circ f=id_V$ y $f\circ f^{-1}=id_W$.
#### Cambio de bases
>[!note] Cambio de bases
>Sean $V$ y $W$ dos $K$-espacios vectoriales de dimensión finita. Sean $B_1,B_1'$ bases de $V$ y $B_2,B_2'$ bases de $W$. Entonces
>$$|f|_{B_1'B_2'}=C(B_2,B_2')\cdot|f|_{B_1B_2}\cdot C(B_1',B_1)$$
##### Demostración
Se tiene que $f=id_W\circ f\circ id_V$. Aplicando la multiplicación de matrices de TLs para obtener la matriz de la composición de TLs, y considerando que la matriz de la TL identidad en un par de bases coincide con la matriz de cambio de base entre las mismas, obtenemos
$$|f|_{B_1'B_2'}=|id_W\circ f\circ id_V|_{B_1'B_2'}=|id_W\circ f|_{B_1B_2'}|id_V|_{B_1'B_1}$$
$$=|id_W|_{B_2B_2'}\cdot|f|_{B_1B_2}\cdot|id_V|_{B_1'B_1}=C(B_2,B_2')\cdot|f|_{B_1B_2}\cdot C(B_1',B_1)$$
que es lo que se quería probar.
### Función inversa de una transformación lineal biyectiva (isomorfismo)
>[!note] Función inversa de un isomorfismo
>Sean $V$ y $W$ dos $K$-espacios vectoriales y sea $f:V\to W$ una TL. Si $f$ es un isomorfismo, entonces $f^{-1}:W\to V$ es una TL (que resulta ser un isomorfismo).
#### Demostración
Sean $w,w'\in W$. Como $f$ es un isomorfismo, existen únicos $v,v'\in V$ tales que $w=f(v)$ y $w'=f(v')$. Entonces
$$f^{-1}(w+w')=f^{-1}(f(v)+f(v'))=f^{-1}(f(v+v'))$$
$$=v+v'=f^{-1}(w)+f^{-1}(w')$$
Dados $w\in W$ y $\lambda\in K$, existe un único $v\in V$ tal que $w=f(v)$. Entonces
$$f^{-1}(\lambda\cdot w)=f^{-1}(\lambda\cdot f(v))=f^{-1}(f(\lambda\cdot v))=\lambda\cdot v=\lambda(f^{-1}(w))$$
Luego, $f^{-1}$ es una TL y es claro que es biyectiva.
## Rango de una matriz
Conociendo la relación entre matrices y transformaciones lineales, podemos estudiar un nuevo invariante asociado a una matriz: su **rango**.
### Rango columna
Sean $V$ y $W$ dos $K$-espacios vectoriales tales que $\text{dim}V=,$ y $\text{dim}W=n$, y sean $B_1$ y $B_2$ bases de $V$ y $W$ respectivamente. Sea $f:V\to W$ una TL. Consideremos la matriz de $f$ en las bases $B_1$ y $B_2$ dada por sus columnas:
$$|f|_{B_1B_2}=(C_1|...|C_m)\in K^{n\times m}$$
Si $B_1=\{v_1,...,v_m\}$, entonces $\text{Im}(f)=<f(v_1),...,f(v_m)>$. Tomando coordenadas en la base $B_2$ se obtiene un subespacio $T\subset K^n$ dado por $T=<(f(v_1))_{B_2},...,(f(v_m))_{B_2}>=<C_1,...,C_m>$. Como tomar coordenadas en una base es un isomorfismo, se tiene que
$$\text{dim}(\text{Im}(f))=\text{dim}<C_1,...,C_m>$$
Entonces, definimos
>[!note] Rango columna
>Sea $A\in K^{n\times m}$. Se llama **rango columna** de $A$, y se nota $\text{rg}_C(A)$, a la dimensión del subespacio de $K^n$ generado por las columnas de $A$, es decir, si $A=(C_1|...|C_m)$, entonces $rg_C(A)=\text{dim}<C_1,...,C_m>$.
#### Importancia
Mediante el cálculo del rango columna de una matriz $A$ es posible obtener la dimensión del subespacio de soluciones del sistema lineal homogéneo asociado a $A$.
#### Observación: rango y nulidad
Sea $A\in K^{n\times m}$ y sea $S=\{x\in K^m: A\cdot x=0\}$. Entonces $\text{dim}S=m-\text{rg}_C(A)$.
En efecto, consideramos la TL asociada a $A$, $f_A:K^m\to K^n$ definida por $f_A(x)=A\cdot x$. Entonces $A=|f_A|_{EE'}$ (donde $E$ y $E'$ son las bases canónicas de $K^m$ y $K^n$ respectivamente) y $S=\text{Nu}(f_A)$. Entonces
$$\text{dim}S=\text{dim}(\text{Nu}(f_A))=m-\text{rg}_C(A)$$
#### Ejemplo
Sea $A\in\mathbb{R}^{3\times3},A=\begin{pmatrix} 1 & -2 & 3 \\ -1 & 2 & 1 \\ 1 & -2 & 4\end{pmatrix}$, y sea $S0\{x\in\mathbb{R}^3:A\cdot x=0\}$. Entonces:
$$\text{dim}S=3-\text{rg}_C(A)=3-\text{dim}<(1,-1,1),(-2,2,-2),(3,1,4)>=3-2=1$$
### Rango fila
Considerando el subespacio generado por las filas de una matriz en lugar del generado por las columnas, puede darse la siguiente definición análoga a la de rango columna.
>[!note] Rango fila
>Sea $A\in K^{n\times m}$. Se define el **rango fila** de $A$, y se nota $\text{rg}_F(A)$, como la dimensión del subespacio de $K^m$ generado por las filas de $A$. Es decir, si $A=\begin{pmatrix} F_1 \\ \vdots \\ F_n \end{pmatrix}$, entonces $\text{rg}_F(A)=\text{dim}<F_1,...,F_n>$.
#### Observación
Sea $A\in K^{n\times m}$. Entonces $\text{rg}_F(A)=\text{rg}_C(A^t)$.
### Lema de Invarianza del Rango
>[!note] Lema de Invarianza del Rango
>Sea $A\in K^{n\times m}$. Sean $C\in\text{GL}(n,K)$ y $D\in\text{GL}(m,K)$. Entonces
>$$\text{rg}_C(A)=\text{rg}_C(C\cdot A\cdot D).$$
#### Demostración
Sea $f_A:K^m\to K^n$ la TL inducida por la multiplicación a izquierda por $A$. Si $E$ y $E'$ son las bases canónicas de $K^m$ y $K^n$ respectivamente, se tiene que $|f_A|_{EE'}=A$ y por lo tanto, $\text{rg}_C(A)=\text{dim}(\text{Im}(f_A))$.
Puesto que $D\in\text{GL}(m,K)$, existe una base $B_1$ de $K^m$ tal que $D=C(B_1,E)$, y como $C\in\text{GL}(n,K)$, existe una base $B_2$ de $K^n$ tal que $C=C(E',B_2)$. Entonces
$$C\cdot A\cdot D=C(E',B_2)\cdot|f_A|_{EE'}\cdot C(B_1,E)=|f_a|_{B_1B_2}$$
de donde $\text{rg}_C(C\cdot A\cdot D)=\text{dim}(\text{Im}(f_A))=\text{rg}_C(A)$.
### Lema de transformación de matrices a la Forma Normal de Jordan
Multiplicando a $A$ por matrices inversibles convenientes se puede obtener una matriz tal que su rango y el de su transpuesta son fáciles de comparar.

>[!note] Lema de transformación de matrices a la Forma Normal de Jordan
>Sea $A\in K^{n\times m}-\{0\}$. Entonces existen $k\in\mathbb{N},1\leq k\leq\text{min}\{n,m\}$, y sus matrices $C\in\text{GL}(n,K)$ y $D\in\text{GL}(m,K)$ tales que
>$$(C\cdot A\cdot D)_{ij}=\begin{cases} 0 & \text{si }i\neq j \\ 1 & \text{si }i=j\leq k \\ 0 & \text{si }i=j>k\end{cases}$$

^379a11

Este lema establece que **cualquier matriz puede ser transformada en una matriz de Jordan** (una matriz que tiene 1s en la diagonal y 0s en todas las demás posiciones) mediante la multiplicación por matrices invertibles a la izquierda y a la derecha. En tu lema, C y D son matrices invertibles, y A es la matriz original. Por lo tanto, la matriz A puede ser transformada en una matriz de Jordan, que es precisamente lo que el Lema de la Forma Normal de Jordan establece.
>[!warning] Versión simplificada de la Forma Normal de Jordan
>Cabe destacar que este lema es una versión simplificada del concepto de la forma normal de Jordan, que generalmente se aplica a matrices cuadradas y puede incluir valores distintos de 1 en la diagonal.
#### Demostración
Consideramos la TL $f_A:K^m\to K^n$ inducida por la multiplicación a izquierda por $A$. Sea $\{v_1,...,v_s\}$ una base de $\text{Nu}(f_A)$ y sean $w_1,...,w_{m-s}\in K^m$ tales que
$$B_1=\{w_1,...,w_{m-s},v_1,...,v_s\}$$
es una base de $K^m$ (si $\text{Nu}(f_A)=\{0\},s=0$ se toma una base $B_1$ cualquiera de $K^m$).
Entonces $\{f_A(w_1),...,f_A(w_{m-s})\}$ es una base de $\text{Im}(f_A)$ y puede extenderse a una base de $K^n$. Sean $z_1,...,z_{n-m+s}\in K^n$ tales que
$$B_2=\{f_A(w_1),...,f_A(w_{m-s}),z_1,...,z_{n-m+s}\}$$
es una base de $K^n$.
Se tiene que
$$(|f_A|_{B_1B_2})_{ij}=\begin{cases}0 & \text{si }i\neq j \\ 1 & \text{si }i=j\leq m-s \\ 0 & \text{si }i=j>m-s\end{cases}$$
Observamos que
$$|f_A|_{B_1B_2}=C(E',B_2)\cdot|f_A|_{EE'}\cdot C(B_1,E)=C\cdot A\cdot D$$
donde $C=C(E',B_2)\in\text{GL}(n,K)$ y $D=C(B_1,E)\in\text{GL}(m,K)$.
### Propiedad de igualdad de rangos columna-fila
>[!note] Rango columna es igual al Rango fila
>Sea $A\in K^{n\times m}$. Entonces $\text{rg}_C(A)=\text{rg}_F(A)$
#### Demostración
Como es claro el resultado para $A=0$. Dada $A\in K^{n\times m}-\{0\}$, por el [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#^379a11|lema de transformación de matrices]], existen matrices $C\in\text{GL}(n,K),D\in\text{GL}(m,K)$ y $k\in\mathbb{N}$, tales que:
$$(C\cdot A\cdot D)_{ij}=\begin{cases} 0 & \text{si }i\neq j \\ 1 & \text{si }i=j\leq k \\ 0 & \text{si }i=j>k\end{cases}$$
Por el [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#Lema de Invarianza del Rango|lema de invarianza del rango]] se tiene que $\text{rg}_C(A)=\text{rg}_C(C\cdot A\cdot D)$, y es claro que $\text{rg}_C(C\cdot A\cdot D)=k$.
Por otro lado, transponiendo obtenemos
$$(C^t\cdot A^t\cdot D^t)_{ij}=\begin{cases} 0 & \text{si }i\neq j \\ 1 & \text{si }i=j\leq k \\ 0 & \text{si }i=j>k\end{cases}$$
con $D^t\in\text{GL}(m,K)$ y $C^t\in\text{GL}(n,K)$, de donde $\text{rg}_C(A^t)=\text{rg}_C(D^t\cdot A^t\cdot C^t)=k$.
En consecuencia: 
$$\text{rg}_F(A)=\text{rg}_C(A^t)=\text{rg}_C(D^t\cdot A^t\cdot C^t)=k=\text{rg}_C(A)$$
### Rango de la matriz
>[!note] Rango de la matriz
>Sea $A\in K^{n\times m}$. Al número $\text{rg}_C(A)=\text{rg}_F(A)$ lo llamaremos el **rango de la matriz** $A$, y lo notaremos $\text{rg}(A)$.
#### Propiedad de dimensión de soluciones
>[!note] Dimensión de soluciones
>Sea $A\in K^{n\times m}$ y sea $S=\{x\in K^m:A\cdot x=0\}$. Entonces $\text{dim}S=m-\text{rg}(A)$.

Esto significa que ==la dimensión del espacio de soluciones de un sistema lineal homogéneo es igual a la cantidad de incógnitas menos la cantidad de ecuaciones independientes==.
## Equivalencia de matrices
>[!note] Matrices equivalentes
>Sean $A,B\in K^{n\times m}$. Se dice que $A$ es equivalente a $B$, y se nota $A\equiv B$, si existen matrices $C\in\text{GL}(n,K)$ y $D\in\text{GL}(m,K)$ tales que $A=C\cdot B\cdot D$.
### El rango y la equivalencia (verificar en los libros)
>[!note] 
>Sean $A,B\in K^{n\times m}$. Entonces $A\equiv B\iff\text{rg}(A)=\text{rg}(B)$.
#### Demostración
##### $\implies$
Es consecuencia del [[Teórica 2 - Transformaciones lineales. Aritmética de punto flotante. Número de condición#Lema de Invarianza del Rango|lema de invarianza de rango]].
##### $\Longleftarrow$
Supongamos que $\text{rg}(A)=\text{rg}(B)=k$. Entonces existen matrices $C_1,C_2\in\text{GL}(n,K)$ y $D_1,D_2\in\text{GL}(m,K)$ tales que
$$(C_1\cdot A\cdot D_1)_{ij}=\begin{cases}0 & \text{si }i\neq j \\ 1 & \text{si }i=j\leq k \\ 0 & \text{si }i=j>k \end{cases}\quad\land\quad(C_2\cdot B\cdot D_2)_{ij}=\begin{cases}0 & \text{si }i\neq j \\ 1 & \text{si }i=j\leq k \\ 0 & \text{si }i=j>k\end{cases}$$
En consecuencia, $C_1\cdot A\cdot D_1=C_2\cdot B\cdot D_2$, de donde
$$A=(C_1^{-1}\cdot C_2)\cdot B\cdot(D_2\cdot D_1^{-1})=C\cdot B\cdot D$$
con $C=C_1^{-1}\cdot C_2\in\text{GL}(n,K)$ y $D=D_2\cdot D_1^{-1}\in\text{GL}(m,K)$.
Luego $A\equiv B$.
### Equivalencia matricial mediante Transformaciones Lineales
Esta propiedad caracteriza matrices equivalentes mediante transformaciones lineales: dos matrices son equivalentes si y sólo si son las matrices de una misma transformación lineal en distintas bases.
>[!note] Equivalencia matricial mediante transformaciones lineales
>Sean $A,B\in K^{n\times m}$. Entonces $A\equiv B$ *si y sólo si* existe una transformación lineal $f:K^m\to K^n$ y bases $B_1,B_1'$ de $K^m$ y $B_2,B_2'$ de $K^m$ tales que $|f|_{B_1B_2}=A$ y $|f|_{B_1'B_2'}=B$.
#### Demostración
##### $\Longleftarrow$
Se deduce de la proposición anterior, considerando que $\text{rg}(A)=\text{dim}(\text{Im}(f))=\text{rg}(B)$.
##### $\implies$
Si consideramos la TL $f:K^m\to K^n$ definida por $f(x)=B\cdot x$. Entonces $B=|f|_{EE'}$, donde $E$ y $E'$ son las bases canónicas de $K^m$ y $K^n$ respectivamente.
Por definición, si $A\equiv B$, existen matrices inversibles $C$ y $D$ tales que $A=C\cdot B\cdot D$. Sea $B_1$ base de $K^m$ tal que $D=C(B_1,E)$ y sea $B_2$ base de $K^n$ tal que $C=(E',B_2)$. Entonces
$$A=C\cdot B\cdot D=C(E',B_2)|f|_{EE'}C(B_1,E)=|f|_{B_1B_2}$$
## Espacios vectoriales de transformaciones lineales
Fijados dos $K$-espacios vectoriales $V$ y $W$, tiene sentido considerar el conjunto de todas las TLs de $V$ en $W$. **Es importante entender la estructura de estos conjuntos de transformaciones lineales**.
### Definición de Espacios de Transformaciones Lineales
>[!note] Espacios de Transformaciones Lineales
>Sean $V$ y $W$ dos $K$-espacios vectoriales. Definimos
>$$\text{Hom}_K(V,W)=\{f:V\to W\text{ tal que }f\ \text{es una transformación lineal}\}$$

Es necesario definir una operación en $\text{Hom}_K(V,W)$ y una acción de $K$ en $\text{Hom}_K(V,W)$ que lo convierten en un $K$-espacio vectorial.
#### Suma
>[!note] Suma
>Dadas $f,g\in\text{Hom}_K(V,W)$, se define $f+g$ como
>$$(f+g)(x)=f(x)+g(x)\qquad\forall x\in V$$
##### Demostración de grupo abeliano
Viendo que $f+g\in\text{Hom}_K(V,W)$, con lo cual $+$ resulta una operación en $\text{Hom}_K(V,W)$
- Es claro que $f+g:V\to W$
- $f+g$ es una TL, pues para cada $x,y\in V$ se tiene que $$(f+g)(x+y)=f(x+y)+g(x+y)=f(x)+f(y)+g(x)+g(y)$$$$=(f(x)+g(x))+(f(y)+g(y))=(f+g)(x)+(f+g)(y)$$Por otro lado, para cada $\mu\in K$ y cada $x\in V$ vale $$(f+g)(\mu\cdot x)=f(\mu\cdot x)+g(\mu\cdot x)=\mu\cdot f(x)+\mu\cdot g(x)$$$$=\mu\cdot(f(x)+g(x))=\mu\cdot(f+g)(x)$$
Es facil verificar que $(\text{Hom}_K(V,W),+)$ es un grupo abeliano.
#### Producto por escalares
>[!note] Producto por escalares
>Dados $f\in\text{Hom}_K(V,W)$ y $\lambda\in K$ se define $(\lambda\cdot f):V\to W$ como
>$$(\lambda\cdot f)(x)=\lambda\cdot f(x)\qquad \forall x\in V$$
##### Demostración
Veamos que $\lambda\cdot f\in\text{Hom}_K(V,W)$, y por lo tanto, $\cdot$ es una acción de $K$ en $\text{Hom}_K(V,W)$:
- Por definición, $\lambda\cdot f:V\to W$.
- $\lambda\cdot f$ es una transformación lineal, pues
	- para todo par de elementos $x,y\in V$: $$(\lambda\cdot f)(x+y)=\lambda\cdot(f(x+y))=\lambda\cdot(f(x)+f(y))=\lambda\cdot f(x)+\lambda\cdot f(y)$$$$=(\lambda\cdot f)(x)+(\lambda\cdot f)(y)$$
	- Para todo $\mu\in K$ y todo $x\in V$: $$(\lambda\cdot f)(\mu\cdot x)=\lambda\cdot(f(\mu\cdot x))=\lambda\cdot(\mu\cdot f(x))=(\lambda\cdot\mu)\cdot f(x)$$$$=\mu\cdot(\lambda\cdot f(x))=\mu\cdot((\lambda\cdot f)(x))$$
#### Propiedades
Además se cumplen las siguientes propiedades: Si $\lambda,\mu\in K$ y $f,g\in\text{Hom}_K(V,W)$
1. $\lambda\cdot(f+g)=\lambda\cdot f+\lambda\cdot g$
2. $(\lambda+\mu)\cdot f=\lambda\cdot f+\mu\cdot f$
3. $1\cdot f=f$
4. $(\lambda\cdot\mu)\cdot f=\lambda\cdot(\mu\cdot f)$
#### Isomorfismo entre el Espacio de Homomorfismos y el Espacio de Matrices
>[!note] Isomorfismo entre el Espacio de Homomorfismos y el Espacio de Matrices
>$(\text{Hom}_K(V,W),+,\cdot)$ es un $K$-espacio vectorial.
>En el caso que ambos $V$ y $W$ son $K$-espacios vectoriales de dimensión finita, $\text{dim}V=n$ y $\text{dim}W=m$, este $K$-espacio vectorial resulta ser isomorfo $K^{m\times n}$
#### Isomorfismo entre TLs y Espacio de Matrices respecto a Bases fijas
>[!note] Isomorfismo entre TLs y Espacio de Matrices respecto a Bases fijas
>Sean $V$ y $W$ dos $K$-espacios vectoriales de dimensión finita, con $\text{dim}V=n$ y $\text{dim}W=m$. Sean $B$ y $B'$ bases de $V$ y $W$ respectivamente. Entonces la función $T:\text{Hom}_K(V,W)\to K^{m\times n}$ definida por $T(f)=|f|_{BB'}$ es un isomorfismo.
##### Demostración
Supongamos que $B=\{v_1,...,v_n\}$ y $B'=\{w_1,...,w_m\}$.
###### $T$ es una transformación lineal
Sean $f,g\in\text{Hom}_K(V,W)$. Por definición, $T(f+g)=|f+g|_{BB'}$. Observemos que la $j$-ésima columna de esta matriz es
$$((f+g)(v_j))_{B'}=(f(v_j)+g(v_j))_{B'}=(f(v_j))_{B'}+(g(v_j))_{B'}$$
es decir, es la suma de las $j$-ésimas columnas de $|f|_{BB'}$ y $|g|_{BB'}$.
Luego, $|f+g|_{BB'}=|f|_{BB'}+|g|_{BB'}$ o, equivalentemente, $T(f+g)=T(f)+T(g)$.
En forma análoga se prueba que $T(\lambda\cdot f)=\lambda\cdot T(f)$.
###### $T$ es un isomorfismo
- **T es monomorfismo**: Seas $F\in\text{Hom}_K(V,W)$ tal que $T(f)=0$, es decir, $|f|_{BB'}=0$. Entonces, $\text{Im}(f)=\{0\}$, de donde $f\equiv0$.
- **T es epimorfismo**: Sea $A\in K^{m\times n}$. Consideramos $f_A:V\to W$ definida por $(f_A(x))_{B'}=(A\cdot(x)_B^t)^t$ para cada $x\in V$.
De este modo, se tiene que $f_A\in\text{Hom}_K(V,W)$ y $T(f_A)=|f_A|_{BB'}=A$.
# Números de Máquina
## Introducción
En Álgebra Lineal Computacional es importante en los algoritmos su **confiabilidad** y su **costo computacional**.
>[!warning] Importancia de la confiabilidad y el costo computacional de los algoritmos
>Puesto que en la práctica no es posible trabajar con precisión infinita, es imprescindible estudiar los **errores numéricos** y su propagación.
>A su vez, puesto que las operaciones realizadas en el algoritmo conllevan una complejidad de cálculo computacional, es importante estudiar el impacto en el tiempo de ejecución.

Para entender los errores numéricos primero debemos definir que son los números máquina.
>[!note] Definición informal de números máquina
>Números que la máquina puede representar de *forma exacta*.

A partir de esta definición, conocemos que los números máquina son finitos, de modo que la mayoría de los números reales no se pueden obtener sin cometer errores.
Para un número real positivo, representados en la base 10 (habitual), este es
$$x=\sum_{-\infty}^{k}\tilde{b}_i10^i, 0\leq\tilde{b}_i\leq9,\tilde{b}_k\neq0$$
La escritura de este modo es única, excepto en desarrollos periódicos del número 9. Por ejemplo $0.999...=1$.
## Bases generales
Todo $x\in\mathbb{R},x\neq0$ puede representarse en una base $b\in\mathbb{N},b\geq2$, de la forma
$$(x)_b=\text{sg}(x)\tilde{b}_k\tilde{b}_{k-1}...\tilde{b}_0\tilde{b}_{-1}\tilde{b}_{-2}...$$
donde los dígitos verifican $0\leq\tilde{b}_i\leq b-1,\text{sg}(x)$ el signo de $x$ y donde por convención numeramos los índices para que el punto -o coma- se ubique entre $\tilde{b}_0$ y $\tilde{b}_{-1}$.
Si obviamos el signo, la representación numérica se obtiene de la escritura de $x$ como
$$x=\tilde{b}_kb^k+\tilde{b}_{k-1}b^{k-1}+...+\tilde{b}_0b^0+\tilde{b}_{-1}b^{-1}+\tilde{b}_{-2}b^{-2}+...$$
que es única si descartamos los desarrollos infinitos de la forma $\tilde{b}_j=b-1$ para todo $j\leq l$ con $l\in\mathbb{Z}$ fijo.
El uso del punto *fijo* determina la división entre la *parte entera* $[x]\in\mathbb{Z}$ y la diferencia $0\leq x-[x]<1$. En efecto, por un lado, para todo $m\in\mathbb{N}$
$$\sum_0^m\tilde{b}_ib^i\in\mathbb{Z}$$
y
$$0\leq\sum_{-\infty}^{-1}\tilde{b}_ib^i<\sum_{-\infty}^{-1}\tilde{b}_i(b-1)=(b-1)\frac{1}{b-1}=1$$
La expresión $0.b_0b_1b_2...$ se llama **mantisa** y $e$ **exponente**.
## Números máquina
>[!note] Números máquina
>Números que pueden representarse de forma exacta con las limitaciones presentes en una máquina, dedicando $m\geq1$ de dígitos para la mantisa y $E\geq1$ para el exponente.
### Cantidad de números máquina
Considerando que $b_0\neq0$ tenemos $b-1$ posibilidades para la elección del $b_0$, y $b$ posibilidades para cada uno de los restantes $b_i$, $1\leq i\leq m-1$. La cantidad diferente de mantisas es entonces $(b-1)b^{m-1}$. Análogamente hay $b^E$ distintos exponentes. Considerando los respectivos signos de la mantisa y el exponente, habrá $(b-1)b^{m+E-1}$ diferentes números máquina.
El mayor exponente para una máquina dada es el número
$$E_{\text{max}}=\sum_{j=0}^{E-1}(b-1)b^j=\frac{b^E-1}{b-1}(b-1)=b^E-1$$
y la mayor mantisa
$$M_\text{max}=\sum_{j=1}^m(b-1)b^{-j}=\frac{(b-1)}{b}\sum_{j=0}^{m-1}b^{-j}=\frac{(b-1)}{b}\frac{(b^{-m}-1)^b}{1-b}=1-b^{-m}$$
lo que da
$$M_\text{max}^{E_\text{max}}=(1-b^{-m})b^{b^E-1}$$
como mayor número de máquina (análogamente se obtiene el menor como $-M_\text{max}^{E_\text{max}}$).
Cualquier mayor producirá un desborde **overflow**. El menor número positivo de máquina es obviamente $b^{-E_\text{max}-1}$. Todo número $x,0<x<b^{-E_\text{max}-1}$ genera un desborde por debajo **underflow**.
### Escalas de los números de máquina
Las máquinas trabajan en base 2 y se utilizan 64 bits para almacenar un número. De estos 64 bits disponibles, dos se reservan para los signos de la mantisa y el exponente.
Considerando que $2^{10}=1024\approx10^3$:
$$M_\text{max}^{E_\text{max}}=(1-2^{-52})2^{2^{10}-1}\approx2^{1000}\approx(2^{10})^{100}\approx10^{300}$$
En la escala pequeña puede trabajar (haciendo una cuenta similar para $b^{-E_\text{max}-1}$) con números del orden $10^{-300}$. Estas cantidades son razonables para describir la mayoría de las magnitudes que aparecen en las disciplinas científicas.
>[!note] Rango de valores de máquina
>Las máquinas de 64-bits pueden trabajar con valores absolutos desde $10^{-300}$ hasta $10^{300}$
#### Redondeo y truncamiento de valores
==Para aquellos números $x\in\mathbb{R}$ que no son de máquina, se elige el número máquina más cercano para representarlo (*redondeo*) o el inmediato inferior (*truncado*).== Dado un real $x$ denominamos $\text{fl}(x)$ al número de máquina elegido para representarlo. Se utiliza $\text{fl}$ para indicar que trabajamos en *punto flotante*.
>[!tip] Movimiento del punto
>El "mover" el punto de acuerdo a la base y el valor del primer decimal, se realiza con la finalidad de mantener uniforme el error relativo a la hora de almacenar $x$ a lo largo de todas las escalas.
### Errores Relativos y Absolutos
>[!note] Error absoluto
>En una magnitud escalar $x_v$, el error absoluto se define como
>$$e_\text{abs}=|x_a-x_v|$$
>donde $x_a$ es el valor aproximado de $x_v$.

En general (y para $x_v\neq0$) estaremos interesados en el error relativo.
>[!note] Error relativo
>El error relativo se define como
>$$e_\text{rel}=\frac{|x_a-x_v|}{|x_v|}$$
>O sea, el tamaño del error absoluto respecto del tamaño del valor exacto $x_v$.
>Para magnitudes vectoriales o matriciales reemplazamos el módulo por una norma.
### Distribución de los números de máquina
Los números máquina no presentan una distribución uniforme. Para describir dos números de máquina consecutivos $x_1<x_2$ asumiendo, para simplificar, que $0<x_1$ y que se escribe:
$$(x_1)_b=0.b_0b_1b_2,,,b_{m-1}b^e$$
Si analizamos esto, el próximo número de máquina (salvo que haya overflow) será
$$(x_2)_b=0-b_0b_1b_2...(b_{m-1}+1)b^e$$
salvo que $b_{m-1}=b-1$ en cuyo caso
$$(x_2)_b=0.b_0b_1b_2...(b_{m-2}+1)0b^e$$
salvo que también $b_{m-2}=b-1$ en cuyo caso habrá que repetir el argumento y eventualmente -en caso de que $b_i=b-1$ para todo $0\leq i\leq m-1$
$$(x_2)_b=0.1b^{e+1}$$
En cualquier caso observamos que
$$x_2-x_1=b^{-m}b^e$$
*no es constante* al variar el exponente $e$, pero sí lo es para cada exponente fijo dado.
>[!tip] Consecuencia en la distribución de los números máquina
>Si $x\in\mathbb{R},\exists x_1<x_2$ números de máquina tales que $x_1\leq x\leq x_2$ (o sea, $x$ está dentro de las escalas manejadas por la máquina) entonces
>$$|x-\text{fl}_\text{redond}(x)|\leq\frac{1}{2}b^{e-m}$$
>$$|x-\text{fl}_\text{trunc}(x)|\leq b^{e-m}$$
>para los casos de redondeo y truncado respectivamente.
### Error relativo y precisión de la máquina
Analizando el caso de redondeo, y por lo tanto considerando que $\text{fl}=\text{fl}_\text{redond}$. Queremos acotar el error relativo al almacenar un número $x\neq0$. Observando que
$$b^{e-1}=|0.1b^e|\leq|x|$$
resulta
$$\left|\frac{x-\text{fl}(x)}{x}\right|\leq\frac{1}{2}b^{1-m}$$
Notar que el error relativo depende del tamaño de la mantisa. En el caso de $b=2$ y $m=52$ mencionado antes resulta:
$$\left|\frac{x-\text{fl}(x)}{x}\right|\leq2^{-52}\approx10^{-16}$$
este número, propio de cada máquina, se denomina **precisión de máquina** o $\epsilon$ de la máquina.
#### Ejemplo: Suma
En general $x\neq\text{fl}(x)$. Sin embargo, podemos garantizar que la diferencia verifica $x-\text{fl}(x)=\mu_xx$ con $\frac{|x-\text{fl}(x)}{|x|}=|\mu_x|\leq\epsilon$.
Si analizamos qué pasa con estos errores durante la suma, distinguiendo esta operación realizada por la máquina como $\oplus$.
Para analizar la suma $x+y$, primero debemos definir las siguientes operaciones:
$$x\to\text{fl}(x)\qquad y\to\text{fl}(y)$$
luego suma $\text{fl}(x)+\text{fl}(y)$ y finalmente almacena el resultado:
$$\text{fl}(x)+\text{fl}(y)\to\text{fl}(\text{fl}(x)+\text{fl}(y))$$
Analizamos el comportamiento del error resultante:
$$\text{fl}(x)=x(1+\mu_x),\qquad\text{fl}(y)=y(1+\mu_y)$$
Por lo tanto,
$$\text{fl}(\text{fl}(x)+\text{fl}(y))=(\text{fl}(x)+\text{fl}(y))(1+\mu_z)$$
y en definitiva
$$x\oplus y=\text{fl}(\text{fl}(x)+\text{fl}(y))=(x(1+\mu_x)+y(1+\mu_y))(1+\mu_z)$$
Vemos entonces que si $0<x,y$ es posible acotar
$$|x\oplus y-(x+y)|\leq(x+y)2\mu+O(\epsilon^2)$$
con $|\mu|\leq\epsilon$. O sea, preservamos de algún modo el tamaño del error relativo. Eso no es posible si restamos números similares.
##### Ejemplo
Sean $b=10,m=4$, la precisión es $\epsilon=\frac{1}{2}10^{-3}=0.0005$. Si $x=125,49$ e $y=125,31$ tenemos respectivamente $x-y=0.18\ominus y=0.2$ por lo que el error relativo es:
$$\left|\frac{x-y-x\ominus y}{x-y}\right|=\frac{0.02}{0.18}\approx0.11$$
O sea, a pesar de que $\epsilon\approx10^{-3}$, el error en la cuenta anterior es del $11\%$. En conclusión, en una simple operación podemos perder dígitos significativos.
>[!warning] Eliminación catastrófica
>Es una regla general **evitar restas de números similares**.

>[!tip] MitigacionCancelacionCatastrofica
>**Cuándo usarla**: Cuando una función o fórmula matemática involucra la resta de dos números casi idénticos de igual signo ($x \approx y \implies x - y \approx 0$).
>**Estrategias resolutivas**:
>1. **Racionalización**: En expresiones con raíces $\sqrt{a} - \sqrt{b}$, multiplicar y dividir por el conjugado $(\sqrt{a} + \sqrt{b})$ transformando la resta en suma.
>2. **Relaciones de Cardano-Vieta en cuadráticas**: Si $ax^2 + bx + c = 0$, calcular la raíz segura sumando números del mismo signo $x_1 = \frac{-b - \operatorname{sign}(b)\sqrt{b^2-4ac}}{2a}$ y despejar la conflictiva mediante $x_2 = \frac{c}{a x_1}$.
>3. **Desarrollo de Taylor**: Para $x \approx 0$ en funciones como $1 - \cos(x)$ o $e^x - 1$, usar su aproximación polinomial centrada en 0.
>**Error típico**: Evaluar directamente la resta confiando en la doble precisión de la computadora; se pueden perder todos los dígitos significativos en un solo paso.
## Condición y estabilidad
### Condicionalidad
Un **problema bien condicionado** es aquel que reacciona benignamente con los errores en los datos. O sea, sus soluciones no varían demasiado al no variar demasiado los datos.
>[!note] Condición de un problema
>Requisitos o restricciones que deben cumplirse para resolver el problema. Estas condiciones pueden ser de varios tipos, como limitaciones de tiempo, recursos disponibles, reglas específicas que deben seguirse, entre otros.

Es un concepto intrínseco al problema, más allá del algoritmo de resolución. Es esencial su condicionamiento para obtener una resolución precisa mediante un algoritmo, siempre que este no incremente desproporcionadamente los errores inherentes a los datos.
>[!tip] La condición es un criterio de resolución
>La condición establecida para un problema funciona como un criterio, de modo que los algoritmos que lo resuelven son clasificados en **estables** e **inestables** según cumplan o no con este criterio. Por lo tanto, la estabilidad es algo *intrínseco del algoritmo*.
### Número de condición
Si consideramos el problema de evaluar una función en el valor $x_0\neq0$. Si por alguna razón (error de medición o redondeo) modificamos el dato a evaluar $x_0$ en una cierta magnitud pequeña $h$ entonces el error relativo que cometemos es (asumiendo que la función es continuamente diferenciable)
$$\frac{f(x_0+h)-f(x_0)}{f(x_0)}=\frac{hf'(\eta)}{f(x_0)}$$
para cierto $\eta$ intermedio entre $x_0$ y $x_0+h$. Esto indica que para $h\sim0$
$$\frac{hf'(\eta)}{f(x_0)}\sim\frac{x_0f'(x_0)}{f(x_0)}\frac{h}{x_0}$$
el error relativo en los datos $\frac{h}{x_0}$ se magnifica en nuestro problema en un factor
$$\frac{|x_0||f'(x_0)|}{|f(x_0)|}$$
llamado **número de condición** de evaluar $f$ en $x_0$. Siguiendo un razonamiento similar vemos que para la versión $f:\mathbb{R}^n\to\mathbb{R}^m$ puede definirse el número de condición de este problema.
#### Ejemplo 1: mal condicionamiento
Un ejemplo concreto y elemental de mal condicionamiento es, como hemos visto, la resta de números similares. Si escribimos $f:\mathbb{R}^2\to\mathbb{R},f(x,y)=x-y$, se tiene $Df(x_0,y_0)=(1,-1),||(1,-1)||_\infty=2,||(x_0,y_0)||_\infty=\text{max}\{|x_0|,|y_0|\}$, luego $\frac{||x_0||||Df(x_0)}{||f(x_0)||}\sim\frac{1}{||x_0-y_0||}$.
#### Ejemplo 2
Si queremos evaluar $\text{tg}(x)$ en un $x_0<\frac{\pi}{2},x_0\sim\frac{\pi}{2}$ vemos que el número de condición
$$\frac{x_0}{\sin{(x_0)}\cos{(x_0)}}$$
cumple que
$$\lim_{x_0\to\frac{\pi}{2}^-}\frac{x_0}{\sin{(x_0)}\cos{(x_0)}}=+\infty$$
En particular, si elegimos $x_0$ tal que $\frac{x_0}{\sin{(x_0)}\cos{(x_0)}}\sim10^{16}$ no esperamos tener ningún dígito significativo en nuestro cálculo de $\text{tg}(x_0)$.
#### Observación: la inestabilidad en iteraciones
La inestabilidad es muy frecuente, fácilmente notable en los resultados obtenidos porque difieren notoriamente de lo esperado. Los casos más extremos aparecen en procesos que deben iterarse en los cuales **el error se acumula exponencialmente** en vez de aditivamente.
#### Ejemplo 3
Supongamos que queremos evaluar $e^{-12}$. Usando lo explicado más arriba, podemos ver que se trata de un problema bien condicionado. Si utilizamos un algoritmo basado en la serie convergente
$$e^{-x}=1-x+\frac{1}{2!}x^2-\frac{1}{3!}x^3...$$
para una evaluación directa en $x=-12$ obtenemos sumando los primeros 50 términos de la serie:
$$e^{-12}\sim6.144189436702122\cdot10^{-6}$$
Si por otro lado modificamos el algoritmo calculado primero $e^{12} sumando los primeros 50 términos de la serie
$$e^x=1+x+\frac{1}{2!}x^2+\frac{1}{3!}x^3...$$
y luego invirtiendo $\frac{1}{e^{12}}$ obtenemos
$$e^{-12}\sim6.144212353328213\cdot10^{-6}$$
Entre estas dos respuestas, la más confiable es la segunda puesto que en el cómputo de la serie alternada, muchos términos "grandes" deben cancelarse mágicamente para producir el resultado "pequeño" $e^{-12}$. Los errores relativos pequeños son proporcionalmente grandes en términos del resultado final.
# Normas vectoriales y sucesiones
## Norma vectorial
Sea el siguiente sistema de $n$ ecuaciones con $n$ incógnitas
$$Ax=b$$
con $A\in\mathbb{R}^{n\times n}$, $x\in\mathbb{R}^n$ y $b\in\mathbb{R}^n$ y nos preguntamos cuánto afectará a la solución un error en el dato $b$. Para poder medir el error, primero definimos como medimos los vectores de $\mathbb{R}^n$. Una forma posible es mediante la norma del vector, o sea,
$$||x||_2=\sqrt{x_1^2+\cdots+x_n^2}$$
También podemos trabajar con otras medidas. Por ejemplo, podemos decir que un vector es "chico" si lo son todas sus componentes y tomar entonces como medida de $x$ la siguiente, llamada **norma infinito**,
$$||x||_\infty=\max_{1\leq i\leq n}|x_i|$$
Otra elección natural es la **norma uno**,
$$||x||_1=|x_1|+\cdots+|x_n|$$
o más en general la **norma p**,
$$||x||_p=(|x_1|^p+\cdots+|x_n|^p)^\frac{1}{p}$$
con $1\leq p<\infty$.
Todas estas formas de medir resultan equivalentes en el sentido de que, si $x$ es "chico" en una de las normas entonces lo es en cualquier otra, puesto que una norma es mayor que la otra salvo una constante que depende sólo de $n$. Por ejemplo, utilizando la [[Clase 1 - Ecuaciones e Inecuaciones. Distancia. Producto Interno. Rectas y Planos#^1e5c5f|desigualdad de Schwartz]] se obtiene
$$||x||_1\leq\sqrt{n}||x||_2$$
y por otra parte, es fácil ver que,
$$||x||_2\leq||x||_1$$
También se verifica fácilmente que
$$||x||_\infty\leq||x||_1\leq n||x||_\infty$$

>[!note] Norma de un vector
>Una norma en $\mathbb{R}^n$ es una manera de asignar a cada $x$ un número $||x||$ de tal forma que se verifiquen las siguientes propiedades, análogas a las que cumple la longitud usual,
>1. $||x||\geq0,\forall x\in\mathbb{R}^n$
>2. $||x||=0\iff x=0$
>3. $||\lambda x||=|\lambda|\ ||x||,\forall x\in\mathbb{R},\forall x\in\mathbb{R}^n$
>4. $||x+y||\leq||x||+||y||,\forall x\in\mathbb{R}^n,\forall y\in\mathbb{R}^n$ (desigualdad triangular).
## Sucesiones
Para medir el error y para analizar la convergencia de una sucesión elegimos la norma que nos resulte más conveniente en cada caso. Esto está justificado por el hecho de que todas las normas son equivalentes.
>[!tip] Equivalencia de norma y convergencia
>La convergencia en una norma implica la convergencia en cualquier otra.
### Teorema de equivalencia de normas
>[!note] Teorema de equivalencia de normas
>Dadas dos normas en $\mathbb{R}^n$, $||\ ||$ y $||\ ||'$, existen constantes $C_1$ y $C_2$ que dependen sólo de $n$ y de las normas consideradas (en particular, son independientes de $x$) tales que
>$$C_1||x||\leq||x||'\leq C_2||x||\qquad\forall x\in\mathbb{R}^n$$
#### Demostración
Basta ver que una norma cualquiera $||\ ||$ es equivalente a la norma euclídea usual, $||\ ||_2$. Sea $\{e_i\}$ la base canónica de $\mathbb{R}^n$ y definamos la constante $C=\left(\sum_{i=1}^n||e_i||^2\right)^{1/2}$, la cual depende sólo de $n$ y de $||\ ||$. Utilizando las propiedades de la norma y la [[Clase 1 - Ecuaciones e Inecuaciones. Distancia. Producto Interno. Rectas y Planos#^1e5c5f|desigualdad de Schwartz]] obtenemos
$$||x||=\left|\left|\sum_{i=1}^nx_ie_i\right|\right|\leq\sum_{i=1}^n|x_i|||e_i||\leq\left(\sum_{i=1}^n|x_i|^2\right)^{1/2}\left(\sum_{i=1}^n||e_i||^2\right)^{1/2}=C||x||_2$$
Queremos ver ahora que existe una constante $K$ tal que
$$||x||_2\leq K||x||\qquad\forall x\in\mathbb{R}^n$$
Supongamos que un tal $K$ no existe y veamos que se llega a una contradicción. En efecto, si esta desigualdad no se cumple para ningún $K$, en particular, que dado $m\in\mathbb{N}$, existe $y_m\in\mathbb{R}^n$ tal que
$$||y_m||_2\geq2||y_m||$$
y llamado $x_m=y_m/||y_m||_2$ obtenemos $||x_m||_2=1$ y
$$||x_m||\leq\frac{1}{m}$$
pero, toda sucesión acotada en la norma euclídea tiene una subsucesión convergente. Entonces, existe una subsucesión de $(x_m),(x_m')$ tal que
$$||x_m'-x||_2\to0$$
para cierto $x\in\mathbb{R}^n$. Pero entonces por lo obtenido mediante la desigualdad de Schwartz, también vale que
$$||x_m'-x||\to0$$
Por otra parte, por $||y_m||_2\geq2||y_m||$ tenemos que $||x_m'||\to0$ y en consecuencia, por unicidad del límite, resulta $x=0$. Pero observemos finalmente que, por desigualdad triangular,
$$\left|||x_m'||_2-||x||_2\right|\leq||x_m'-x||_2$$
y entonces se llega a la contradicción $1=||x_m'||_2\to||x||_2=0$.
### Norma matricial asociada
>[!note] Norma matricial asociada
>Definimos para una matriz arbitraria $A\in\mathbb{R}^{n\times n}$ y una norma de vectores $||\ ||$ cualquiera, la norma matricial asociada es
>$$||A||=\max_{0\neq x\in\mathbb{R}^n}\frac{||Ax||}{||x||}$$

O sea, la norma de $A$ nos da lo máximo que "agranda" el multiplicar $A$ medido en la norma de vectores dada. Es fácil ver que
$$||A||=\max_{||x||=1}||Ax||$$
y en particular, *esto muestra que el máximo existe, o sea que la norma está bien definida* ($||Ax||$ es un función continua de $x$ y por lo tanto, alcanza su máximo en el conjunto de vectores de norma igual a uno pues éste es cerrado y acotado).
> [!tip] Desigualdad de Norma Matricial Asociada
> De la definición se desprende la siguiente desigualdad que usaremos frecuentemente,
$$||Ax||\leq||A||\cdot||x||,\quad\forall x\in\mathbb{R}^n$$
valiendo la igualdad para algún $x$.

Por otra parte, puede verificarse que $||A||$ es la menor entre todas las constantes $C$ para las cuales vale la desigualdad
$$||Ax||\leq C||x||,\quad\forall x\in\mathbb{R}^n$$
siendo esta otra forma usual de definir la norma matricial.
### Desigualdad de Normas Matriciales
Sean $A,B\in\mathbb K^{n\times n}$ tenemos que:
$$||AB||\leq||A||\cdot||B||$$
## Factor de amplificación del error relativo
Si se reemplaza el dato $b$ por $b+\Delta b$, de la solución $x$ del sistema será modificada de tal forma que tendremos
$$A(x+\Delta x)=(b+\Delta b)$$
o equivalentemente,
$$A\Delta x=\Delta b$$
y nos preguntamos que relación hay entre
$$\frac{||\Delta x||}{||x||}\qquad\text{y}\qquad\frac{||\Delta b||}{||b||}$$
Veamos primero el siguiente ejemplo simple,
$$\begin{pmatrix}
	4,1 & 2,8\\
	9,7 & 6,6
\end{pmatrix}
\begin{pmatrix}
	x_1 \\
	x_2
\end{pmatrix}
=
\begin{pmatrix}
	4,1 \\
	9,7
\end{pmatrix}
$$
La solución es
$$x=\begin{pmatrix} 1 \\ 0 \end{pmatrix}$$
Observemos que
$$||b||_1=13,8\qquad||x||_1=1$$
Si modificamos $b$ poniendo
$$b'=b+\Delta b=\begin{pmatrix} 4,11 \\ 9,70 \end{pmatrix}$$
entonces la solución es
$$x'=x+\Delta x=\begin{pmatrix}0,34 \\ 0,97\end{pmatrix}$$
y se obtiene en consecuencia
$$||\Delta b||_1=0,01\qquad||\Delta x||_1=1,63$$
con lo que el error relativo se amplificó mucho, en efecto,
$$\frac{||\Delta b||_1}{||b||_1}=0,7246\times10^{-3}$$
mientras que
$$\frac{||\Delta x||_1}{||x||_1}=1,63$$
La importancia es entender a qué se debe este comportamiento y poder predecir, dada una matriz $A$, cuál será el **factor de amplificación del error relativo** o, al menos, dar una cota de éste en términos de $A$.
Analicemos primero el caso de una matriz diagonal.
$$\begin{pmatrix}\lambda_1 & 0 \\ 0 & \lambda_2\end{pmatrix}
\begin{pmatrix}x_1 \\ x_2\end{pmatrix}=
\begin{pmatrix} b_1 \\ b_2 \end{pmatrix}$$
La solución es (si $\lambda_1,\lambda_2\neq0)$
$$x_1=\frac{b_1}{\lambda_1}\qquad x_2=\frac{b_2}{\lambda_2}$$
Por ejemplo, si
$$A=\begin{pmatrix} 1000 & 0 \\ 0 & \frac{1}{100}\end{pmatrix}$$
entonces,
$$x_1=\frac{b_1}{1000}\qquad x_2=100b_2$$
Si ponemos $b'=b+\Delta b$ con
$$\Delta b=\begin{pmatrix}0\\\Delta b_2\end{pmatrix}\qquad
b=\begin{pmatrix}b_1 \\ 0\end{pmatrix}$$
entonces,
$$x_1=\frac{b_1}{1000}\qquad \Delta x_2=100\Delta b_2$$
obteniéndose
$$10^5\frac{||\Delta b||_2}{||b||_2}=\frac{||\Delta x||_2}{||x||_2}$$
es decir que el error relativo se multiplicó por $10^5$.
Si en cambio elegimos
$$\Delta b=\begin{pmatrix}\Delta b_1 \\ 0\end{pmatrix}\qquad
b=\begin{pmatrix}0 \\ b_2\end{pmatrix}$$
tenemos entonces,
$$\Delta x_1=\frac{1}{1000}\Delta b_1\qquad x_2=100b_2$$
y en consecuencia,
$$\frac{1}{10^5}\frac{||\Delta b||_2}{||b||_2}=\frac{||\Delta x||_2}{||x||_2}$$
o sea que en este caso el error relativo se redujo en un factor $10^5$. En general, para una matriz diagonal
$$A=\begin{pmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{pmatrix}$$
con $|\lambda_1|>|\lambda_2|$, el error relativo puede multiplicarse por
$$\frac{|\lambda_1|}{|\lambda_2|}$$
en el peor de los casos y por
$$\frac{|\lambda_2|}{|\lambda_1|}$$
en el mejor de los casos.
En general, el error tendrá componentes en cualquier dirección por lo que es de esperar que si $\frac{|\lambda_1|}{|\lambda_2|}$ es grande los errores relativos se amplifiquen.
El mismo análisis puede hacerse en $\mathbb{R}n$. Si $A$ es una matriz diagonal
$$A=\begin{pmatrix}\lambda_1 & \cdots & 0 \\ & \cdots & \\ 0 & \cdots & \lambda_N\end{pmatrix}$$
el error relativo se puede amplificar, a lo sumo, por un factor
$$\frac{|\lambda_\max|}{|\lambda_\min|}$$
siendo $\lambda_\max$ y $\lambda_\min$ los de máximo y mínimo valor absoluto entre los $\lambda_j$ (puesto que $A$ es inversible, $\lambda_\min\neq0$). Este cociente se llama **número de condición** o condicionamiento de $A$ en la norma $||\ ||_2$ y lo denotaremos $\text{Cond}_2(A)$.
## Número de condición
Sea la matriz $A$ simétrica, o sea $a_{ij}=a_{ji}$. En este caso $A$ se puede diagonalizar, es decir, existe una base de autovectores $\{v_1,...,v_n\}$. Además, por ser $A$ simétrica podemos considerar que la base es ortonormal. Entonces, si $Ax=b,A(x+\Delta x)=b+\Delta b$ y
$$x=\sum_{i=1}^n\alpha_iv_i\qquad\Delta x=\sum_{i=1}^n\beta_iv_i$$
tenemos,
$$||x||_2^2=\sum_{i=1}^n\alpha_i^2\qquad||\Delta x||_2^2=\sum_{i=1}^n\beta_i^2$$
y además, si llamamos $\lambda_i$ al autovalor correspondiente a $v_i$,
$$b=\sum_{i=1}^n\alpha_i\lambda_iv_i\qquad\Delta b=\sum_{i=1}^n\beta_i\lambda_iv_i$$
y en consecuencia,
$$||b||_2^2=\sum_{i=1}^n\alpha_i^2\lambda_i^2\qquad||\Delta b||_2^2=\sum_{i=1}^n\beta_i^2\lambda_i^2$$
Entonces, si $\lambda_\max$ y $\lambda_\min$ son los valores de máximo y mínimo valor absoluto respectivamente, obtenemos
$$\frac{||\Delta x||_2^2}{||x||_2^2}=\frac{\sum_{i=1}^n\beta_i^2}{\sum_{i=1}^n\alpha_i^2}\leq\frac{1/|\lambda_\min|^2}{1/|\lambda_\max|^2}\frac{||\Delta b||_2^2}{||b||_2^2}$$
o sea
$$\frac{||\Delta x||_2}{||x||_2}\leq\frac{|\lambda_\max|}{|\lambda_\min|}\frac{||\Delta b||_2}{||b||_2}$$

>[!note] Número de condición
El número $\text{Cond}_2(A)=\frac{|\lambda_\max|}{|\lambda_\min|}$ es una cota para el factor de amplificación del error relativo. Esta cota es la mejor posible pues la desigualdad se convierte en una igualdad para cierta elección de $b$ y $\Delta b$ ($b$ en la dirección correspondiente al máximo autovalor y $\Delta b$ en la correspondiente al mínimo).

Podemos generalizar el número de condición a cualquier matriz si observamos que, siendo $A$ simétrica, la dirección correspondiente al autovalor de máximo valor absoluto nos da la dirección de "máxima expansión", o sea, si miramos el cociente entre la longitud de $Ax$ y la de $x$
$$\frac{||Ax||_2}{||x||_2}$$
éste será máximo entre todos los $x$ cuando $x$ está en la dirección correspondiente a $\lambda_\max$. En efecto, si escribimos $x$ en la base de autovectores $\{v_1,...,v_n\}$
$$x=\sum_{i=1}^n\alpha_iv_i$$
entonces,
$$Ax=\sum_{i=1}^n\alpha_i\lambda_iv_i$$
y de acá resulta que
$$||Ax||_2\leq|\lambda_\max|\cdot||x||_2$$
y tomando $x=v_j$ con $\lambda_j=\lambda_\max$ se verifica la igualdad.
Análogamente, la dirección de "mínima expansión" corresponde a la asociada a $\lambda_\min$, la cual corresponde también a la de "máxima expansión" de la inversa de $A$.
>[!info] Factor de amplificación del error relativo en matrices simétricas
>Está relacionado con los máximos factores de expansión de $A$ y de su inversa.
### Definición del número de condición de matriz
>[!note] Número de condición de matriz
>Sea $A\in\mathbb{R}^{n\times n}$ una matriz inversible y sea $||\cdot||$ una norma en $\mathbb{R}^n$ definimos el número de condición de $A$ como
>$$\text{Cond}(A)=||A||\cdot||A^{-1}||$$

Es claro que $\text{Cond}(A)$ depende de la norma de vectores elegida.
>[!tip] Motivación ([[Práctica 2 - Aritmética de punto flotante. Número de condición#Ejercicio 19|Ejercicio 19]])
>Sean $Ax=b$, tenemos métodos que nos aproximan la solución. Estos métodos nos dan un sistema similar $A\tilde x=\tilde b$.
>Para conocer que tan cerca está el sistema $Ax=b\sim A\tilde x=\tilde b$ es mediante el número de condición.
### Propiedades del número de condición
- $\text{Cond}(A)=\text{Cond}(A^{-1})$
- $\text{Cond}(A)\geq1,\quad\forall A\in\mathbb{R}^{n\times n}$
Para verificar esta segunda propiedad, utilizamos la [[Clase 1 - Ecuaciones e Inecuaciones. Distancia. Producto Interno. Rectas y Planos#^1e5c5f|desigualdad de Cauchy-Schwarz]] y obtenemos
$$1=||I||=||A\cdot A^{-1}||\leq||A||\cdot||A^{-1}||=\text{Cond}(A)$$
### Teorema de Condición de la Matriz y Proporción de Cambio
>[!note] Teorema de Condición de la Matriz y Proporción de Cambio
>Si $A\in\mathbb{R}^{n\times n}$ es inversible, $b$ ,$\Delta b\in\mathbb{R}^n$, $Ax=b$ y $A(x+\Delta x)=b+\Delta b$ entonces,
>$$\frac{||\Delta x||}{||x||}\leq\text{Cond}(A)\frac{||\Delta b||}{||b||}$$
>valiendo la igualdad para alguna elección de $b$ y $\Delta b$.
>Además,
>$$\frac{1}{\text{Cond}(A)}\frac{||\Delta b||}{||b||}\leq\frac{||\Delta x||}{||x||}$$ 
>y nuevamente, vale la igualdad para ciertos $b$ y $\Delta b$.
#### Demostración
Se tiene que $A(\Delta x)=\Delta b$ y entonces,
$$\frac{||\Delta x||}{||x||}\leq\frac{||A^{-1}||\,||\Delta b||}{||x||}=\frac{||A^{-1}||\,||\Delta b||}{||b||}\frac{||b||}{||x||}\leq\frac{||A^{-1}||\,||\Delta b||}{||b||}||A||$$
donde para la última desigualdad hemos usado que $||b||=||Ax||\leq||A||\,||x||$. Por lo tanto, el teorema anterior queda verificado.
Ahora, si elegimos $\Delta b$ tal que $||\Delta x||=||A^{-1}\Delta b||=||A^{-1}||\,||\Delta b||$, $x$ tal que $||Ax||=||A||\,||x||$ (lo que siempre se puede por la definición de la norma matricial) y $b=Ax$ resulta la igualdad en el teorema de Condición de la Matriz.
Ahora, para ver la segunda desigualdad del Teorema de Condición de Matriz, observemos que ésta puede escribirse como:
$$\frac{||\Delta b||}{||b||}\leq\text{Cond}(A)\frac{||\Delta x||}{||x||}$$
la cual, teniendo en cuenta que $\text{Cond}(A)=\text{Cond}(A^{-1})$, es exactamente la primer desigualdad aplicada a $A^{-1}$ con lo que el teorema está demostrado.
### Teorema Condición de Matriz y la Perturbación
>[!note] Teorema Condición de Matriz y la Perturbación
>Si $A\in\mathbb{R}^{n\times n}$ es inversible, $E\in\mathbb{R}^{n\times n}$, $b\in\mathbb{R}^n$, $Ax=b$ y $(A+E)(x+\Delta x)=b$ entonces, llamando $\tilde{x}=x+\Delta x$ tenemos
>$$\frac{||\Delta x||}{||\tilde{x}||}\leq\text{Cond}(A)\frac{||E||}{||A||}$$
>y si
>$$\text{Cond}(A)\frac{||E||}{||A||}\leq\delta<1$$
>entonces
>$$\frac{||\Delta x||}{||x||}\leq\frac{1}{1-\delta}\text{Cond}(A)\frac{||E||}{||A||}$$
#### Demostración
Teniendo que $Ax=b$ y $A\tilde{x}=b-E\tilde{x}$, entonces $-E\tilde{x})A\Delta x$. Por lo tanto,
$$||\Delta x||\leq||A^{-1}||\,||E||\,||\tilde{x}||\leq\text{Cond}(A)\frac{||E||}{||A||}||\tilde{x}||$$
lo que prueba el último teorema.
Ahora, siendo que
$$\frac{||\tilde{x}||}{||x||}=\frac{||x+\Delta x||}{||x||}\leq\frac{||x||+||\Delta x||}{||x||}=1+\frac{||\Delta x||}{||x||}$$
lo cual, junto con la desigualdad anterior implica que
$$\frac{||\Delta x||}{||x||}\leq\text{Cond}(A)\frac{||E||}{||A||}\left(1+\frac{||\Delta x||}{||x||}\right)\leq\text{Cond}(A)\frac{||E||}{||A||}+\delta\frac{||\Delta x||}{||x||}$$
### Radio espectral de una matriz
>[!note] Radio espectral de una matriz
>Dada $A\in\mathbb{R}^{n\times n}$ se llama **radio espectral** de $A$ a
>$$\rho(A)=|\lambda_\max|$$
>siendo $\lambda_\max$ el de máximo valor absoluto entre todos los autovalores de $A$, incluyendo los complejos.
### Norma de Frobenius
Si $A$ es simétrica entonces
$$||A||_2=\rho(A)$$
En general, para $A\in\mathbb{R}^{n\times n}$ arbitraria se tiene
$$||A||_2=\sqrt{\rho(A^tA)}$$
Como $A^tA$ es simétrica, existe una base ortonormal de autovectores $\{v_j\}$. Llamando $\mu_j$ a los autovalores correspondientes tenemos
$$A^tAv_j=\mu_jv_j\qquad j=1,...,n$$
y si $x=\sum_{i=1}^n\alpha_iv_i$ entonces, por la ortonormalidad de los $v_j$ resulta $||x||_2^2=\sum_{i=1}^n\alpha_j^2$ y en consecuencia, considerando que $A^tAx=\sum_{i=1}^n\alpha_j\mu_jv_j$, se tiene que para todo $x\in\mathbb{R}^n$
$$\frac{||Ax||_2^2}{||x||_2^2}=\frac{x^tA^tAx}{||x||_2^2}=\frac{\sum_{i=1}^n\alpha_j^2\mu_j}{\sum_{i=1}^n\alpha_j^2}\leq\rho(A^tA)$$
es decir que
$$||A||_2\leq\sqrt{\rho(A^tA)}$$
y tomando $x=v_j$ con $\mu_j=\mu_\max$ se ve que vale la igualdad.
>[!tip] La excepcional norma de Frobenius
>Las normas no subordinadas de matrices no suelen cumplir con la desigualdad de Cauchy-Schwarz. Pero, la norma de Frobenius es una excepción. Dados $A,B\in\mathbb{K}^{n\times n}$, consideremos el elemento $i,j$ de la matriz producto.
>$$|[AB]_{i,j}|^2\leq\left|\left(\sum_{l=1}^nA_{i,l}B_{l,j}\right)\right|^2\leq\left(\sum_{l=1}^nA_{i,l}^2\right)\left(\sum_{l=1}^nB_{l,j}^2\right)$$
>donde usando la desigualdad de Cauchy-Schwarz en la última desigualdad. Sumando en $i$ y en $j$ se tiene el resultado deseado
>$$||AB||_F^2=\sum_{1\leq i,j\leq n}|[AB]_{i,j}|^2\leq\sum_{1\leq i\leq n}\sum_{1\leq j\leq n}\left(\sum_{l=1}^nA_{i,l}^2\right)\left(\sum_{l=1}^nB_{l,j}^2\right)=||A||_F^2||B||_F^2$$
>esto es
>$$||AB||_F\leq||A||_F||B||_F$$
### Elección de norma matricial y cálculo

| Norma | Nombre | Fórmula de Cálculo Práctico | Tipo | Complejidad |
| :--- | :--- | :--- | :--- | :--- |
| **$\|A\|_1$** | Suma por columnas | $\max_{1 \le j \le n} \sum_{i=1}^m |a_{ij}|$ (máxima suma de valores absolutos de una columna) | Inducida | $O(mn)$ |
| **$\|A\|_\infty$** | Suma por filas | $\max_{1 \le i \le m} \sum_{j=1}^n |a_{ij}|$ (máxima suma de valores absolutos de una fila) | Inducida | $O(mn)$ |
| **$\|A\|_2$** | Espectral | $\sigma_{\max}(A) = \sqrt{\lambda_{\max}(A^T A)}$ | Inducida | $O(n^3)$ |
| **$\|A\|_F$** | Frobenius | $\sqrt{\sum_{i=1}^m \sum_{j=1}^n a_{ij}^2} = \sqrt{\operatorname{tr}(A^T A)}$ | No inducida ($\|I\|=\sqrt{n}$) | $O(mn)$ |

>[!tip] CalculoNormasMatriciales
>**Cuándo usarla**: En ejercicios de parcial para acotar o calcular rápidamente la norma de una matriz sin calcular autovalores ni valores singulares.
>**Pasos**:
>1. Para $\|A\|_1$: sumar el valor absoluto de cada **columna** y tomar el máximo.
>2. Para $\|A\|_\infty$: sumar el valor absoluto de cada **fila** y tomar el máximo.
>3. Si la matriz es diagonal $D$: $\|D\|_1 = \|D\|_\infty = \|D\|_2 = \max_i |d_{ii}|$.
>**Error típico**: Confundir filas con columnas entre norma 1 y norma $\infty$, o calcular $\|A\|_2$ sumando al cuadrado las entradas (eso es Frobenius $\|A\|_F$, no norma 2).

Por ejemplo, para la norma $\infty$ se tiene:
$$||A||_\infty=\max_{1\leq i\leq m}\sum_{j=1}^n|a_{ij}|$$
Para ver esto observemos primero que, para todo $x\in\mathbb{R}^n$ vale
$$||Ax||_\infty=\max_{1\leq i\leq n}\left|\sum_{j=1}^Na_{ij}x_j\right|\leq||x||_\infty\max_{1\leq i\leq n}\sum_{j=1}^N|a_{ij}|$$
Para ver que vale la otra desigualdad, sea $k$ tal que $\sum_{j=1}^n|a_{kj}|$ es máxima y tomemos
$$x=\begin{pmatrix}
	\text{sg}(a_{k1}) \\
	\text{sg}(a_{k2}) \\
	\cdots \\
	\text{sg}(a_{kn})
\end{pmatrix}$$
donde $\text{sg}(a)=1$ si $a\geq0$ y $\text{sg}(a)=-1$ si $a<0$. Entonces
$$(Ax)_k=\sum_{j=1}^n|a_{kj}|$$
y en particular, $||Ax||_\infty\geq\sum_{j=1}^n|a_{kj}|$ y como $||x||_\infty=1$ obtenemos
$$\frac{||A_x||_\infty}{||x||_\infty}\geq\sum_{j=1}^n|a_{kj}|$$
y concluimos que
$$||A||_\infty\geq\max_{1\leq i\leq n}\sum_{j=1}^n|a_{ij}|$$
De manera similar puede verse, aunque no lo demostraremos, que
$$||A||_1=\max_{1\leq j\leq n}\sum_{i=1}^n|a_{ij}|$$
### Importancia práctica del número de condición de una matriz

>[!tip] CotaErrorPerturbacionCondicion
>**Cuándo usarla**: Para estimar o acotar la pérdida de precisión o el error relativo de la solución $x$ en un sistema $Ax=b$ cuando hay errores de redondeo o mediciones con ruido en $b$ o en $A$.
>**Fórmulas clave**:
>1. Perturbación en $b$ ($A(x+\delta x) = b + \delta b$):
>   $$\frac{\|\delta x\|}{\|x\|} \le \kappa(A) \frac{\|\delta b\|}{\|b\|}$$
>2. Perturbación general en $A$ y en $b$ (con $\kappa(A)\frac{\|\delta A\|}{\|A\|} < 1$):
>   $$\frac{\|\delta x\|}{\|x\|} \le \frac{\kappa(A)}{1 - \kappa(A)\frac{\|\delta A\|}{\|A\|}} \left( \frac{\|\delta b\|}{\|b\|} + \frac{\|\delta A\|}{\|A\|} \right)$$
>3. Regla de dígitos significativos:
>   $$\text{Dígitos correctos garantizados} \approx (\text{Dígitos máquina}) - \log_{10}(\kappa(A))$$
>**Error típico**: Multiplicar $\sigma_{\max} \cdot \sigma_{\min}$ en lugar de dividir $\kappa_2(A) = \frac{\sigma_{\max}(A)}{\sigma_{\min}(A)}$ para matrices simétricas/diagonales. Recordar que $\kappa(A) \ge 1$ siempre.

**Referencias y Conexiones**:
- *Se usa en*: [[Práctica 2 - Aritmética de punto flotante. Número de condición]] Ej. 13 a 20.
- *Criterio de diseño*: Es la justificación teórica de por qué en [[Teórica 3 - Sistemas lineales y factorización#Descomposición $LU=PA$ (pivoteo parcial)]] es obligatorio hacer pivoteo parcial (para mantener acotado el condicionamiento de los factores triangulares).

Supongamos una $A$ matriz singular y, para cierto $b$, el sistema $Ax=b$ tiene alguna solución, entonces tendrá infinitas y estas formarán una variedad lineal de $\mathbb{R}^n$. O sea, sin cambiar $b$, se pueden obtener soluciones tan distantes como se quiera. En otras palabras, en este caso tendríamos $\Delta b=0$ mientras que $\Delta x$ sería arbitrariamente grande.
#### Ejemplo
Sea $\epsilon\sim0$ entonces la matriz
$$A=\begin{pmatrix}1 & 1+\epsilon \\ 1-\epsilon & 1 \end{pmatrix}$$
está cerca de la matriz singular
$$B=\begin{pmatrix}1 & 1 \\ 1 & 1\end{pmatrix}$$
y en este caso
$$A^{-1}=\epsilon^{-2}\begin{pmatrix}1 & -1-\epsilon \\ -1+\epsilon & 1\end{pmatrix}$$
entonces,
$$||A||_\infty=2+\epsilon\qquad||A^{-1}||_\infty=(2+\epsilon)\epsilon^{-2}$$
y en consecuencia,
$$\text{Cond}_\infty(A)=\left(\frac{2+\epsilon}{\epsilon}\right)^2>\frac{4}{\epsilon^2}$$

>[!warning] Distancia de matrices singulares
>Deben entenderse en forma relativa al tamaño de $A$.

En este ejemplo tenemos no sólo que $||A-B||_\infty$ es chica sino que lo es en relación con $||A||_\infty$. En efecto,
$$\frac{||A-B||_\infty}{||A||_\infty}=\frac{\epsilon}{2+\epsilon}<\frac{\epsilon}{2}$$
#### Ejemplo 2
Que una matriz esté cerca de ser singular no tiene nada que ver con el tamaño del determinante. Para verificar esto, consideremos $\epsilon\sim0$, la matriz
$$A\begin{pmatrix}\epsilon & 0 \\ 0 & \epsilon\end{pmatrix}$$
tiene determinante muy pequeño, pero es una matriz muy buena en cuanto a la propagación de errores relativos pues $\text{Cond}(A)=1$.
En cambio,
$$A=\begin{pmatrix}1 & 0 \\ 0 & \frac{1}{\epsilon}\end{pmatrix}$$
tiene determinante grande pero, en las normas $2$, $1$ o $\infty$, $\text{Cond}(A)=\frac{1}{\epsilon}$.
### Relación del número de condición con su distancia relativa a las matrices singulares
>[!note] Teorema de la Cota Inferior del número de condición de una matriz
>Dadas $A\in\mathbb{R}^{n\times n}$ inversible y una norma de vectores cualquiera se tiene
>$$\frac{1}{\text{Cond}(A)}=\inf_{B\ \text{singular}}\frac{||A-B||}{||A||}$$
#### Demostración
Sea $B\in\mathbb{R}^{n\times n}$ una matriz singular y tomemos $x\neq0$ tal que $Bx=0$. Entonces,
$$||x||=||A^{-1}(A-B)x||\leq||A^{-1}||\,||A-B||\,||x||$$
y en consecuencia
$$1\leq||A^{-1}||\,||A-B||$$
lo cual muestra que
$$\frac{1}{\text{Cond}(A)}\leq\frac{||A-B||}{||A||}\qquad\forall B\in\mathbb{R}^{n\times n}\ \text{singular}$$
Entonces, para concluir el teorema, falta ver que hay una $B$ singular para la cual vale la igualdad. Para esto, sean $y$ tal que $||A^{-1}y||=||A^{-1}||\,||y||$ y $x$ tal que $Ax=y$. Como $y$ puede tomarse con norma arbitraria, lo elegimos de tal forma que $||y||=\frac{1}{||A^{-1}||}$ y en consecuencia $||x||=1$.
Sea $z$ un vector tal que
$$z^tx=1\qquad\text{y}\qquad z^tu\leq1,\forall u\in\mathbb{R}^n\ \text{tal que}\ ||u||=1$$
La existencia de un $z$ podemos demostrarla al observar que, intuitivamente, si analizamos el significado geométrico: los $u$ que verifican $z^tu=1$ forman un hiperplano (una variedad lineal de dimensión $n-1$). Por lo tanto, que haya un $z$ verificando estas dos ecuaciones, significa que hay un hiperplano que pasa por $x$ y que deja a la bola unitaria $\mathcal{B}_1=\{u\in\mathbb{R}^n:||u||\leq1\}$ de un lado. La existencia de tal hiperplano es clara si se considera que, para toda norma, $\mathcal{B}_1$ es un conjunto convexo y que $x$ está en el borde de éste. Observemos también que, en el caso de la norma $||\ ||_2$, se tiene que $z=x$.
Ahora, definimos $B=A-yz^t$ y vemos que esta matriz es singular y cumple con la igualdad que queríamos. En efecto,
$$Bx=Ax-yz^tx=y-y=0$$
y por lo tanto, $B$ es singular.
Por otra parte, $||A-B||=||yz^t||$, pero por la definición de $z$ tenemos que $|z^tu|\leq1$ para todo $u$ tal que $||u||=1$ puesto que, si $z^tu<0$ entonces, $|z^tu|=-z^tu=z^t(-u)\leq1$ ya que $||-u||=1$. Entonces
$$||yz^tu||=||y||\,|z^tu|\leq\frac{1}{||A^{-1}||}\qquad\forall u\in\mathbb{R}^n\ \text{tal que}\ ||u||=1$$
y por lo tanto,
$$||A-B||\leq\frac{1}{||A^{-1}||}$$
concluyendo con la demostración.