# 🌊 Esquema Numérico para la Ecuación 2D Sine-Gordon Bifraccional Disipativa

## 🎯 Resumen del Proyecto

[cite_start]Este repositorio contiene la implementación computacional en **MATLAB** de un esquema de **diferencias finitas** de segundo orden para resolver la Ecuación 2D Sine-Gordon (FSG) con doble amortiguamiento (disipación) y operadores de **Riesz espacio-fraccionales** de diferente orden ($\alpha, \beta \in (1, 2)$)[cite: 714].

[cite_start]El objetivo central del trabajo fue desarrollar un modelo numérico que preservara la estructura de disipación (o conservación, en el caso sin amortiguamiento) de la energía del sistema, tal como ocurre en el modelo continuo[cite: 82].

**Publicación Asociada:**
[cite_start]El análisis y los resultados de este trabajo han sido publicados en la revista *Fractal and Fractional*[cite: 1, 14, 15].

## ⚛️ Modelo Matemático (Ecuación Sine-Gordon Fraccional 2D)

Se investigó la siguiente Ecuación Diferencial Parcial no lineal:

$$
\frac{\partial^{2}u(x,y,t)}{\partial t^{2}}+\gamma\frac{\partial u(x,y,t)}{\partial t}-\lambda\Delta^{\alpha,\beta}u(x,y,t)=-\phi(x,y)sin~u(x,y,t) +F(x,y,t)-G^{\prime}(u(x,y,t))
$$
[cite_start]donde $\Delta^{\alpha,\beta}u = \frac{\partial^{\alpha}u}{\partial|x|^{\alpha}}+\frac{\partial^{\beta}u}{\partial|y|^{\beta}}$ es el **Operador Laplaciano Fraccional de orden $(\alpha, \beta)$** de Riesz[cite: 120].

## 🛠️ Metodología Numérica

* [cite_start]**Esquema Temporal:** Se empleó un método de **Crank-Nicolson** de dos pasos para la discretización temporal[cite: 86, 722].
* [cite_start]**Aproximación Espacial:** Los operadores fraccionales de Riesz se aproximaron mediante **diferencias centradas de orden fraccional**, lo que ofrece una precisión de segundo orden en espacio y tiempo[cite: 87, 721, 722].
* [cite_start]**Solución del Sistema:** El sistema no lineal resultante se resolvió mediante un esquema de **punto fijo** (*fixed-point iteration*)[cite: 429, 725, 729].

### 🔒 Propiedades Teóricas Demostradas

El esquema numérico satisface rigurosamente las siguientes propiedades, fundamentales para la fiabilidad de la solución:

1.  [cite_start]**Conservación/Disipación de Energía Discreta (Preservación de Estructura):** Se probó que la energía discreta se conserva (cuando $\gamma=0$ y $F \equiv 0$) o disipa (cuando $\gamma>0$) a lo largo del tiempo, replicando el comportamiento del sistema continuo (Teorema 2)[cite: 31, 82, 724].
2.  [cite_start]**Consistencia:** El esquema es de **segundo orden** en el error de truncamiento ($\mathcal{O}(h^2 + \tau^2)$)[cite: 87, 322, 726].
3.  [cite_start]**Estabilidad Condicional y Convergencia:** La estabilidad y la convergencia del esquema son de **segundo orden en la norma $L^2$**[cite: 87, 412, 700, 727].

## 🧪 Resultados de Simulación

Las simulaciones en MATLAB confirmaron la validez del esquema:

* [cite_start]**Verificación:** El error $L^2$ entre la solución numérica y una solución analítica conocida se mantuvo por debajo de $1.6 \times 10^{-2}$[cite: 440].
* [cite_start]**Preservación de Energía:** Se validó la propiedad de conservación/disipación de energía, coincidiendo con los resultados del Teorema 2[cite: 32, 578].
* [cite_start]**Efectos Fraccionales:** Al variar las órdenes fraccionales ($\alpha, \beta$), se observó que la **amplitud de la onda aumenta** significativamente a medida que $\alpha$ y $\beta$ se alejan del caso entero ($\alpha=\beta=2$)[cite: 487].

## 💻 Ejecución y Código

[cite_start]El código de simulación está escrito en **MATLAB (versión 2024)**[cite: 774].

* El archivo principal para la simulación es `simulation.m`.
* El código completo de las funciones auxiliares se encuentra en el Apéndice A del artículo.

---

## 👥 Autores y Contacto

**Autores:**
[cite_start]Dagoberto Mares-Rincón, Siegfried Macías, Jorge E. Macías-Díaz, José A. Guerrero-Díaz-de-León, and Tassos Bountis[cite: 6].

[![Gmail Badge](https://img.shields.io/badge/-dagobertomares0@gmail.com-c14438?style=flat&logo=Gmail&logoColor=white&link=mailto:dagobertomares0@gmail.com)](mailto:dagobertomares0@gmail.com) - 
[![Linkedin Badge](https://img.shields.io/badge/-dagobertomares-0072b1?style=flat&logo=Linkedin&logoColor=white&link=https://www.linkedin.com/in/dagoberto-mares/)](https://www.linkedin.com/in/dagoberto-mares/)
