"En bash o terminal, asegúrate de tener HMMER instalado para ejecutar este script. Puedes instalarlo usando:"
"sudo apt-get install hmmer  # Ubuntu/Debian"
"brew install hmmer          # macOS"

import os
import subprocess
import shutil
import re

# Configuración de rutas de carpetas de entrada y salida
CARPETA_PFAM = "modelos_hmm_pfam"       
CARPETA_FASTA = "uniprot_fastas"        
CARPETA_RESULTADOS = "resultados_hmmer"  


def buscar_binario_wsl(nombre_binario):
    """Busca un binario en el PATH o en las rutas estándar de Ubuntu WSL."""
    ruta = shutil.which(nombre_binario)
    if ruta:
        return ruta

    rutas_comunes = [
        f"/usr/bin/{nombre_binario}",
        f"/usr/local/bin/{nombre_binario}",
        f"/bin/{nombre_binario}"
    ]
    for r in rutas_comunes:
        if os.path.exists(r) and os.access(r, os.X_OK):
            return r
    return None


def inicializar_entorno():
    """Crea las carpetas de salida necesarias."""
    os.makedirs(CARPETA_RESULTADOS, exist_ok=True)


def preparar_base_datos_pfam(carpeta_pfam, ruta_hmmpress):
    """Concatena todos los modelos HMM individuales en una base de datos única y la prensa."""
    db_combinada = os.path.join(carpeta_pfam, "pfam_db.hmm")
    
    archivos_hmm = [
        os.path.join(carpeta_pfam, f) 
        for f in os.listdir(carpeta_pfam) 
        if f.endswith(".hmm") and not f.startswith("pfam_db")
    ]
    
    if not archivos_hmm:
        raise FileNotFoundError(f"No se encontraron archivos .hmm en la carpeta '{carpeta_pfam}'")
    
    print(f"-> Unificando {len(archivos_hmm)} perfiles HMM en una base de datos...")
    with open(db_combinada, "wb") as f_salida:
        for ruta_hmm in archivos_hmm:
            with open(ruta_hmm, "rb") as f_entrada:
                f_salida.write(f_entrada.read())
                
    print("-> Indexando base de datos con hmmpress en WSL...")
    try:
        subprocess.run([ruta_hmmpress, "-f", db_combinada], check=True, stdout=subprocess.DEVNULL)
        print("¡Base de datos de Pfam lista y optimizada!")
        return db_combinada
    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar hmmpress dentro de WSL: {e}")
        return None


def ejecutar_hmmscan(db_pfam, carpeta_fasta, carpeta_salida, ruta_hmmscan):
    """Escanea cada archivo FASTA contra la base de datos de Pfam unificada."""
    archivos_fasta = [
        os.path.join(carpeta_fasta, f) 
        for f in os.listdir(carpeta_fasta) 
        if f.endswith(".fasta") or f.endswith(".fa")
    ]
    
    if not archivos_fasta:
        print(f"No se encontraron archivos FASTA en '{carpeta_fasta}'.")
        return

    print(f"\nIniciando escaneo de {len(archivos_fasta)} archivos de proteínas...")
    for ruta_fasta in archivos_fasta:
        nombre_base = os.path.splitext(os.path.basename(ruta_fasta))[0]
        archivo_reporte = os.path.join(carpeta_salida, f"resultado_{nombre_base}.txt")
        archivo_tabla = os.path.join(carpeta_salida, f"tabla_{nombre_base}.tbl")
        
        print(f"   Escaneando {nombre_base} contra perfiles Pfam...", end="", flush=True)
        comando = [
            ruta_hmmscan,
            "--domtblout", archivo_tabla,
            "-o", archivo_reporte,
            db_pfam,
            ruta_fasta
        ]
        try:
            subprocess.run(comando, check=True)
            print(" OK")
        except subprocess.CalledProcessError as e:
            print(f" Error al procesar {nombre_base}: {e}")


def limpiar_archivos_temporales(carpeta_pfam):
    """Elimina de forma segura la base de datos indexada y los archivos auxiliares binarios."""
    print("\nIniciando limpieza de archivos binarios e índices temporales...")
    extensiones_basura = [".hmm", ".hmm.h3f", ".hmm.h3i", ".hmm.h3m", ".hmm.h3p"]
    archivos_eliminados = 0

    for ext in extensiones_basura:
        ruta_basura = os.path.join(carpeta_pfam, f"pfam_db{ext}")
        if os.path.exists(ruta_basura):
            try:
                os.remove(ruta_basura)
                archivos_eliminados += 1
            except Exception as e:
                print(f"   [Aviso] No se pudo eliminar {ruta_basura}: {e}")

    if archivos_eliminados > 0:
        print(f"¡Limpieza completada! Se eliminaron {archivos_eliminados} archivos residuales.")
    else:
        print("No se encontraron archivos temporales para limpiar.")


if __name__ == "__main__":
    bin_hmmscan = buscar_binario_wsl("hmmscan")
    bin_hmmpress = buscar_binario_wsl("hmmpress")

    if not bin_hmmscan or not bin_hmmpress:
        print("\n[ERROR CRÍTICO]: Faltan herramientas esenciales de HMMER en WSL Ubuntu.")
        print("Instálalas corriendo: sudo apt-get install -y hmmer\n")
    else:
        print("HMMER3 validado de forma exitosa en el entorno WSL Ubuntu.")
        inicializar_entorno()
        try:
            ruta_db_lista = preparar_base_datos_pfam(CARPETA_PFAM, bin_hmmpress)
            if ruta_db_lista:
                ejecutar_hmmscan(ruta_db_lista, CARPETA_FASTA, CARPETA_RESULTADOS, bin_hmmscan)
                limpiar_archivos_temporales(CARPETA_PFAM)
                print("\n[PIPELINE COMPLETADO CON ÉXITO] Procesamiento de HMMER finalizado.")
        except Exception as e:
            print(f"\nSe detuvo el proceso debido a un error inesperado: {e}")
