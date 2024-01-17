#! /bin/bash
# Descripción: Obten las métricas del fichero de logs
# Autor: Maite Bernaus (maite.bernaus@gmail.com)
# Fecha 26/12/2023
# Versión : 1.0


# Verificar la cantidad de parámetros
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 archivo_entrada archivo_salida"
    exit 1
fi

# Obtener nombres de archivos de los parámetros
input_file="$1"
output_file="$2"

# Comprobar si el archivo de salida ya existe
if [ ! -f "$output_file" ]; then
    # Si no existe, agregar el encabezado
    echo -e "Seccion\tTiempo (s)\tMemoria (KB)" >> "$output_file"
fi

# Procesar el archivo con awk
awk '
    /^### / {
        # Guardar el nombre de la sección
        section_name = substr($0, 5)
    }
    /^[0-9]+(\.[0-9]+)?[[:space:]][0-9]+$/ {
        # Extraer tiempo y memoria y escribir en el archivo de salida
        time = $1
        memory = $2
        print section_name "\t" time "\t" memory >> "'$output_file'"
    }
' "$input_file"

echo "Resumen creado en $output_file"

