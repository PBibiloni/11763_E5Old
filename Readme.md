# 11763 Procesamiento de Imagen Médica - Ejercicio 5

Se plantean los siguientes objetivos:

1. Representación de una traslación 3D
2. Representación de una rotación axial
    * Multiplicación de cuaterniones
    * Cálculo de rotaciones axiales mediante su expresión en cuaterniones
3. Corregistro mediante landmarks: supongamos que en dos imágenes, hemos manualmente localizado unos puntos de interés. Buscaremos el ajuste que *mejor* envía los puntos de interés de una a la otra.
    * Calcular los residuos cuadráticos y el error cuadrático medio
    * Encontrar una inicialización
    * Encontrar una buena transformación rígida (¿la mejor?)
    * Estudiar si nuestra representación contiene parámetros redundantes

Pedro Bibiloni
Universitat de les Illes Balears
