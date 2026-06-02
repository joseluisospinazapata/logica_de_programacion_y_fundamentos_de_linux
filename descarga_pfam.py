import os
import time
import gzip
from datetime import datetime
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def configurar_sesion():
    """Configura una sesión de requests con reintentos automáticos."""
    sesion = requests.Session()
    reintentos = Retry(
        total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504]
    )
    sesion.mount("https://", HTTPAdapter(max_retries=reintentos))
    return sesion


def construir_url_dinamica(acc_familia) -> str:
    """Construye de manera dinámica la URL de descarga para cada ID de Pfam."""
    dominio_base = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam"
    endpoint_recurso = "?annotation=hmm"
    return f"{dominio_base}/{acc_familia}{endpoint_recurso}"


def registrar_log(mensaje, archivo_log="registro_descarga.log"):
    """Escribe una línea con marca de tiempo en el archivo de log."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    linea = f"[{timestamp}] {mensaje}\n"
    with open(archivo_log, "a", encoding="utf-8") as f:
        f.write(linea)


def es_archivo_gz(ruta_archivo) -> bool:
    """Verifica si un archivo está comprimido en formato GZIP leyendo sus bytes mágicos."""
    try:
        with open(ruta_archivo, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except Exception:
        return False


def manejar_descompresion(ruta_archivo, acc) -> str:
    """Detecta si el archivo es .gz, lo descomprime y limpia el archivo comprimido."""
    if es_archivo_gz(ruta_archivo):
        ruta_descomprimido = ruta_archivo[:-3]  # Reemplazará el archivo directamente a .hmm plano
        ruta_temporal_gz = ruta_archivo + ".gz"
        
        # Renombrar temporalmente el archivo para operar la descompresión
        os.rename(ruta_archivo, ruta_temporal_gz)
        
        try:
            with gzip.open(ruta_temporal_gz, "rb") as f_entrada:
                contenido_descomprimido = f_entrada.read()
                
            with open(ruta_descomprimido, "wb") as f_salida:
                f_salida.write(contenido_descomprimido)
                
            os.remove(ruta_temporal_gz)
            return " (Archivo .gz detectado y descomprimido exitosamente)"
        except Exception as e:
            # En caso de error, intentar restaurar el archivo original
            if os.path.exists(ruta_temporal_gz) and not os.path.exists(ruta_descomprimido):
                os.rename(ruta_temporal_gz, ruta_descomprimido)
            return f" (Error al descomprimir .gz: {str(e)})"
            
    return " (Texto plano recibido, no requiere descompresión)"


def descargar_modelos_hmm(
    lista_familias, carpeta_salida="modelos_hmm_pfam"
) -> None:
    """Descarga los perfiles HMM, maneja descompresión .gz dinámica y genera logs."""
    os.makedirs(carpeta_salida, exist_ok=True)
    sesion = configurar_sesion()

    mensaje_inicio = f"Iniciando descarga masiva de {len(lista_familias)} familias de proteínas."
    print(mensaje_inicio)
    registrar_log(mensaje_inicio)

    for i, acc in enumerate(lista_familias, 1):
        acc = acc.strip()
        if not acc:
            continue

        ruta_archivo = os.path.join(carpeta_salida, f"{acc}.gz")
        url_dinamica = construir_url_dinamica(acc)

        # Evitar re-descargar si el archivo ya existe
        if os.path.exists(ruta_archivo):
            msg_saltado = f"[{i}/{len(lista_familias)}] Saltado: {acc}.gz ya existe."
            print(msg_saltado)
            registrar_log(msg_saltado)
            continue

        print(f"[{i}/{len(lista_familias)}] Descargando {acc}... ", end="")

        try:
            respuesta = sesion.get(url_dinamica, timeout=20)

            if respuesta.status_code == 200:
                # Guardar el contenido crudo inicialmente
                with open(ruta_archivo, "wb") as f:
                    f.write(respuesta.content)
                
                # Procesar descompresión dinámica si aplica
                resultado_compresion = manejar_descompresion(ruta_archivo, acc)
                
                print(f"OK{resultado_compresion}")
                registrar_log(f"ÉXITO - ID: {acc} | URL: {url_dinamica}{resultado_compresion}")
            else:
                print(f"Error {respuesta.status_code}")
                registrar_log(
                    f"ERROR {respuesta.status_code} - ID: {acc} | URL: {url_dinamica}"
                )

        except Exception as e:
            print("Falló por red")
            registrar_log(
                f"FALLO DE CONEXIÓN - ID: {acc} | Error: {str(e)} | URL: {url_dinamica}"
            )

        # Pausa de cortesía para el servidor de la API
        time.sleep(0.5)

    mensaje_fin = "Proceso de descarga finalizado."
    print(mensaje_fin)
    registrar_log(mensaje_fin)


if __name__ == "__main__":
    # Lista de las 38 familias solicitadas
    lista_familias = [
        "PF00069", "PF07714", "PF00071", "PF00017", "PF00018", "PF00096",
        "PF00046", "PF00013", "PF00170", "PF00447", "PF00076", "PF00270",
        "PF00015", "PF00039", "PF00005", "PF00161", "PF00520", "PF01241",
        "PF00155", "PF00171", "PF00085", "PF00118", "PF00043", "PF00400",
        "PF00023", "PF00020", "PF00092", "PF00012", "PF00183", "PF00182",
        "PF00072", "PF07719", "PF00675", "PF04542", "PF00047", "PF00008",
        "PF00028", "PF00041"
    ]

    descargar_modelos_hmm(lista_familias)
