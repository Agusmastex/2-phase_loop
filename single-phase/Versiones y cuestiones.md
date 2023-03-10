# Features de las versiones


# Transient Temperature

2 - Usa dos for loops para el espacio y tiempo  
3 - Usa matrices de diferenciación  
4 - Guarda todos los estados pasados, T es una matriz n+1 x nt

En 1 y 2 no había nodos en las fronteras, y había un nodo fantasma más allá de cada frontera.  
En 3 y 4 el problema se resuelve para todo el vector T, y después se actualizan las fronteras a la fuerza.  
Todos usan upwinding.

## Notas sobre el comportamiento y la estabilidad
### Método explícito
Cuando la malla es gruesa: la onda se propaga de manera más suave (difusión numérica artificial), y la temperatura máxima es menor  
Cuando la malla es fina: la onda se propaga de manera más abrupta y rectilínea; y se alcanza la temperatura máxima calculada analíticamente.  
Teniendo en cuenta la solución analítica de la ecuación diferencial, creo que la malla fina es más exacta.

### Método implícito
Actualizar la frontera a la fuerza no funciona con el método implícito.  
La frontera de inflow, tiende a cero. Forzar el valor actualizando no cambia nada.  
Cuando la malla es gruesa: se comporta como el método explícito cuando su malla es gruesa.  
Cuando la malla es fina: se genera un frente de onda que pasa por encima de la temperatura máxima y se propaga hasta el final. Tras eso, el sistema vuelve a alcanzar el estado estacionario esperado.

# Simple Transient Velocity

Se eliminó el término de Boussinesq para desacoplar las ecuaciones de energía y movimiento, de modo a estudiar mejor el método de resolución de la velocidad en estado transitorio. La velocidad inicial siempre es cero.

PressureDriven - Se utilizó dp/dx = -ΔP/L constante  
FreeFall &emsp; &emsp;&emsp;- Se utilizó dp/dx = 0  
OnOffPump&emsp;&emsp;- Se utilizó dp/dx como una función escalón unitario

Al hallar la caída de presión con Vz dada para flujo forzado ascendente, esta ΔP es la caída de presión para la cual las velocidades de estado estacionario en flujo forzado ascendente y en caída libre son de igual magnitud. ¿Por qué?

# Transient Averaged

Calcula el perfil de temperaturas transitorio y la evolución de Vz en el tiempo con Δp dada. Utiliza la ecuación de movimiento promediada en el espacio. Cambia la matriz de diferenciación cuando el flujo se invierte. Ejemplos: flow_reversal.gif

To do: Construir una función que tome los parámetros importantes (th,tf,Δp,β) y genere un gif. Loopear a lo largo de los valores de β (o Δp, o th) para encontrar el valor crítico por encima del cual el flujo no se invierte y se mantiene por convección natural. (Obs: este valor está cerca de β = 0.03 para los valores de prueba utilizados)

# ODE Solvers

 - steady: Resuelve los perfiles de presión y temperaturas con Vz dada
 - transient_temperature: Resuelve el perfil de temperaturas transitorio con Vz dada
 - transient: Resuelve los perfiles de presión y temperaturas transitorios con Vz dada
