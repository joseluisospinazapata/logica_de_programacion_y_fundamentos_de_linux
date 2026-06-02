import json
import os
import time
from datetime import datetime
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

# Nombres de archivos y directorios locales
CARPETA_FASTA = "uniprot_fastas"
ARCHIVO_VERSIONES = "control_versiones_uniprot.json"
ARCHIVO_LOG_ERRORES = "errores_uniprot.log"


def configurar_sesion():
    """Configura una sesión de requests con reintentos para evitar microcortes."""
    sesion = requests.Session()
    reintentos = Retry(
        total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504]
    )
    sesion.mount("https://", HTTPAdapter(max_retries=reintentos))
    return sesion


def construir_url_dinamica(lista_ids) -> str:
    """Construye dinámicamente la URL para obtener el FASTA de UniProt."""
    return f"https://rest.uniprot.org/uniprotkb/{lista_ids}.fasta"


def registrar_error(id_uniprot, motivo, url):
    """Escribe en el log de errores si un ID falla o si el FASTA está vacío."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    linea = f"[{timestamp}] ID: {id_uniprot} | MOTIVO: {motivo} | URL: {url}\n"
    with open(ARCHIVO_LOG_ERRORES, "a", encoding="utf-8") as f:
        f.write(linea)


def cargar_control_versiones():
    """Carga el historial de ETags locales para validar cambios futuros."""
    if os.path.exists(ARCHIVO_VERSIONES):
        with open(ARCHIVO_VERSIONES, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def guardar_control_versiones(datos):
    """Guarda el historial actualizado de ETags en el archivo JSON."""
    with open(ARCHIVO_VERSIONES, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=4)


def procesar_uniprot_masivo(lista_ids):
    """Descarga los FASTAs, controla versiones futuras y gestiona el log de errores."""
    os.makedirs(CARPETA_FASTA, exist_ok=True)
    sesion = configurar_sesion()
    historial_versiones = cargar_control_versiones()

    print(f"Iniciando el procesamiento de {len(lista_ids)} IDs de UniProt...")

    for i, id_actual in enumerate(lista_ids, 1):
        id_actual = id_actual.strip()
        if not id_actual:
            continue

        url_dinamica = construir_url_dinamica(id_actual)
        ruta_fasta = os.path.join(CARPETA_FASTA, f"{id_actual}.fasta")

        print(f"[{i}/{len(lista_ids)}] Verificando {id_actual}... ", end="")

        try:
            # Hacemos una petición HEAD primero para comprobar la versión (ETag) sin bajar todo el archivo
            respuesta_head = sesion.head(url_dinamica, timeout=15)

            if respuesta_head.status_code != 200:
                print(f"ERROR HTTP {respuesta_head.status_code}")
                registrar_error(
                    id_actual,
                    f"Código HTTP de error {respuesta_head.status_code}",
                    url_dinamica,
                )
                continue

            # UniProt asigna un ETag único a cada versión de la secuencia
            etag_servidor = respuesta_head.headers.get("ETag", "").replace(
                '"', ""
            )
            etag_local = historial_versiones.get(id_actual)

            # LÓGICA DE VERIFICACIÓN AUTOMATIZADA DE MODIFICACIONES FUTURAS
            if os.path.exists(ruta_fasta) and etag_local:
                if etag_local == etag_servidor:
                    print("Saltado (Secuencia sin cambios en el servidor).")
                    continue
                else:
                    print(
                        "¡ALERTA! La proteína ha sido MODIFICADA en UniProt. Redescargando..."
                    )
            else:
                print("Descargando nueva secuencia... ", end="")

            # Si el archivo no existe o el ETag cambió, descargamos el contenido real
            respuesta_get = sesion.get(url_dinamica, timeout=15)
            contenido_fasta = respuesta_get.text

            # VERIFICACIÓN DE CONTENIDO (Loguea si no contiene la secuencia o está vacío)
            if not contenido_fasta or ">" not in contenido_fasta:
                print("FALLO (Archivo vacío o sin formato FASTA)")
                registrar_error(
                    id_actual,
                    "El archivo descargado no contiene una secuencia FASTA válida",
                    url_dinamica,
                )
                continue

            # Guardar el archivo FASTA localmente
            with open(ruta_fasta, "w", encoding="utf-8") as f_fasta:
                f_fasta.write(contenido_fasta)

            # Actualizar el registro de versión local
            historial_versiones[id_actual] = etag_servidor
            guardar_control_versiones(historial_versiones)
            print("OK")

        except requests.exceptions.RequestException as e:
            print("ERROR DE CONEXIÓN")
            registrar_error(
                id_actual, f"Fallo de red / Timeout: {str(e)}", url_dinamica
            )

        # Pausa de cortesía para la API de UniProt
        time.sleep(0.3)

    print("\nProceso completado.")
    if os.path.exists(ARCHIVO_LOG_ERRORES):
        print(
            f"Nota: Revisa '{ARCHIVO_LOG_ERRORES}' por si algún ID falló la verificación de contenido."
        )


if __name__ == "__main__":
    # Tus 50 identificadores validados de UniProt
    lista_ids = [
        "P00519",
        "P42684",
        "P12931",
        "P06241",
        "P07947",
        "Q06187",
        "P43403",
        "P43405",
        "P62993",
        "P01112",
        "P01116",
        "P01111",
        "P04049",
        "P31749",
        "P28482",
        "P27361",
        "P00533",
        "P21802",
        "P16234",
        "P12956",
        "P29353",
        "P42681",
        "P35222",
        "P62937",
        "P29317",
        "P08047",
        "P15056",
        "P40763",
        "P42224",
        "P15924",
        "P10242",
        "P19838",
        "P11473",
        "P61244",
        "Q9Y2T1",
        "P08107",
        "P0A6Y8",
        "P0A6W5",
        "P0A9Q7",
        "P0A799",
        "P0A7Y4",
        "P0A8V2",
        "P39451",
        "P11142",
        "P13569",
        "P22681",
        "P98160",
        "P12814",
        "Q92793",
        "Q13485",
    ]

    procesar_uniprot_masivo(lista_ids)
