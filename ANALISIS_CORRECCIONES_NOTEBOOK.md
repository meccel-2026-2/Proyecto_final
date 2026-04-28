# Análisis y Correcciones de Cálculos en el Notebook

## Resumen Ejecutivo
He revisado los tres cálculos principales del notebook y encontré **un error crítico** en el cálculo del vector del plano invariante de Laplace. Los otros dos cálculos (energía total y momentum angular) son correctos.

---

## 1. ✅ ENERGÍA TOTAL (CORRECTO)

### Implementación
```python
# Energía cinética
v_2 = np.sum(vs_Raro**2, axis=2)
Ks = 0.5 * ms[:, np.newaxis] * v_2
K = np.sum(Ks, axis=0)

# Energía potencial
U_ij = - G * ms[i] * ms[j] / r_ij_norm
U = np.sum(U_i, axis=0)

# Energía total
E = K + U
```

### Verificación
✓ La energía cinética: $K = \frac{1}{2}\sum_i m_i v_i^2$ es correcta

✓ La energía potencial: $U = -\sum_{i<j} G \frac{m_i m_j}{r_{ij}}$ es correcta

✓ La energía total: $E = K + U$ es correcta

**Conclusión:** Sin problemas. El cálculo es físicamente correcto.

---

## 2. ✅ MOMENTUM ANGULAR TOTAL (CORRECTO)

### Implementación
```python
L_i = np.cross(rs_Raro, P_i)  # donde P_i = m_i * vs_Raro
L_tots_vec = np.sum(L_i, axis=0)
```

### Verificación
✓ Calcula correctamente: $\vec{L}_{tot} = \sum_i m_i \vec{r}_i \times \vec{v}_i$

✓ Como `rs_Raro` contiene posiciones ABSOLUTAS (no relativas), `L_tots_vec` es el momentum angular total absoluto del sistema.

**Conclusión:** Sin problemas. El cálculo es correcto.

---

## 3. ❌ VECTOR DEL PLANO INVARIANTE DE LAPLACE (ERROR CRÍTICO)

### Implementación INCORRECTA (Original)
```python
L_tots_vec = np.sum(L_i, axis=0)  # Momentum angular absoluto

L_tots_sis_vec = L_tots_vec + np.cross(R_cm_vec, P_tots_vec)  # ¡ERROR!
```

### El Problema
El código suma el término del centro de masa **DOS VECES**:

Matemáticamente, el momentum angular total se descompone como:
$$\vec{L}_{tot} = \sum_i m_i \vec{r}_i \times \vec{v}_i = \sum_i m_i \vec{r}_i' \times \vec{v}_i' + M\vec{R}_{CM} \times \vec{V}_{CM}$$

Donde:
- $\vec{r}_i'$ = posición relativa al CM
- $\vec{V}_{CM}$ = velocidad del CM

Puesto que `rs_Raro` son posiciones ABSOLUTAS, `L_tots_vec` ya contiene **ambos términos**:
- El momentum angular relativo al CM: $\sum_i m_i \vec{r}_i' \times \vec{v}_i'$
- El momentum angular del CM: $M\vec{R}_{CM} \times \vec{V}_{CM}$

**Sumar nuevamente el término del CM es incorrecto** porque se cuenta dos veces.

### Implementación CORRECTA
```python
# El plano invariante de Laplace ES simplemente:
L_Laplace = L_tots_vec  # Ya contiene todos los términos correctamente

# Si deseas los componentes por separado:
L_relativo_al_CM = L_tots_vec - np.cross(R_cm_vec, P_tots_vec)
L_del_CM = np.cross(R_cm_vec, P_tots_vec)
```

### Verificación de la Corrección
- $\vec{L}_{Laplace} = \vec{L}_{relativo} + \vec{L}_{CM}$ ✓
- $\vec{L}_{relativo} = \sum_i m_i \vec{r}_i' \times \vec{v}_i'$ ✓  
- $\vec{L}_{CM} = M\vec{R}_{CM} \times \vec{V}_{CM}$ ✓

---

## Correcciones Realizadas

### Cambio 1: Línea 281-286 (Cálculo del momentum angular)
**Antes:** Imprimía solo L_tots_vec sin explicación
**Después:** Se aclara que L_tots_vec es el momentum angular total (plano invariante de Laplace)

### Cambio 2: Línea 330-336 (Descomposición del momentum angular)
**Antes:**
```python
L_tots_sis_vec = L_tots_vec + np.cross(R_cm_vec, P_tots_vec)
L_tots_sis = np.linalg.norm(L_tots_sis_vec, axis=1)
```

**Después:**
```python
# Descomposición correcta
L_relativo_al_CM = L_tots_vec - np.cross(R_cm_vec, P_tots_vec)
L_del_CM = np.cross(R_cm_vec, P_tots_vec)
L_Laplace = L_tots_vec  # Vector del plano invariante

# Impresión de componentes
```

### Cambio 3: Línea 338 (Visualización)
**Antes:** Graficaba `np.linalg.norm(L_tots_sis_vec, axis=1)` sin contexto
**Después:** Grafica la descomposición completa mostrando L_total, L_relativo y L_CM

---

## Implicaciones Físicas

### Antes de la Corrección
- El cálculo del plano invariante de Laplace estaba **incorrecto**
- Los análisis posteriores basados en este vector estarían viciados

### Después de la Corrección  
- El plano invariante de Laplace ahora se calcula correctamente como el momentum angular total
- La conservación de este vector es un principio fundamental en dinámica orbital
- Se puede analizar correctamente la estabilidad orbital de Apophis

---

## Recomendaciones Adicionales

1. **Verificar conservación:** Graficar $\Delta \vec{L} = \vec{L}(t) - \vec{L}(t_0)$ para verificar que sea aproximadamente cero

2. **Análisis virial:** La condición de virial para N-cuerpos es $2\langle K \rangle + \langle U \rangle \approx 0$, que debería verificarse

3. **Documentación:** Agregar comentarios explicando la descomposición del momentum angular en futuras versiones

---

## Referencias Teóricas

- Murray & Dermott - "Solar System Dynamics" - Cap. 2
- Goldstein - "Classical Mechanics" - Cap. 3 (Momento Angular)
- Teoría de N-cuerpos: https://en.wikipedia.org/wiki/N-body_problem
