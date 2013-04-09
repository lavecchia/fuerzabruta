Fuerza Bruta Racional
===================

Requisitos: 
PyMol

Script para búsqueda conformacional de moléculas. El mismo es desarrollado en Python y utiliza librerías de PyMol para la modificación de parámetros moleculares.
Ha sido desarrollado para utilizarlo en sistemas operativos Linux, pero es fácilmente adaptable a otros sistemas compatibles con Python.

Para su utilización, debe crearse en el mismo directorio en donde se encuentra el archivo PDB con la estructura a modificar, un archivo con extensión ".in" con el mismo nombre que este.
En el archivo ".in" se deben especificar en líneas separadas los "atom ID" correspondientes que definen el ángulo diedro del enlace a rotar precedidos por "D" (hace referencia de diedro). Por ejemplo, si se desea rotar el enlace atomo1-atomo2---atomo3-atomo4 se escribe:
"D 1 2 3 4"

Si no se especifica nada luego, los valores que tomara este diedro corresponderán a los mínimos del perfil energético de rotación.
En el caso que se quiera que el diedro tome valores específicos, se debe definir de la siguiente manera:
"D 1 2 3 4 = -130 30 160"
En este ejemplo, este diedro tomará los valores -130º, 30º y 160º.

El script se ejecuta de la siguiente manera:

forzabruta -i estructura.pdb

Este genera, en función de lo definido en el archivo estructura.in, los distintas estructuras (en archivos .pdb) producto de la combinatoria de los valores de cada diedro.
El archivo de salida contendrá el siguiente formato:

estructura_D0_valor0_D1_valor1_..._Dn_valorn.pdb


