"En bash o terminal, asegúrate de tener HMMER instalado para ejecutar este script. Puedes instalarlo usando:"
"sudo apt-get install hmmer  # Ubuntu/Debian"
"brew install hmmer          # macOS"

import os
import subprocess
import shutil

# Configuración de rutas de carpetas de entrada y salida
CARPETA_PFAM = "modelos_hmm_pfam"       # Tu carpeta con los .hmm de Pfam
CARPETA_FASTA = "uniprot_fastas"        # Tu carpeta con las proteínas en FASTA
CARPETA_RESULTADOS = "resultados_hmmer"  # Carpeta donde se guardará todo


def buscar_binario_wsl(nombre_binario):
    """Busca un binario en el PATH o en las rutas estándar de Ubuntu WSL

    para evitar problemas de entorno con VS Code.
    """
    # 1. Intento estándar usando el PATH actual
    ruta = shutil.which(nombre_binario)
    if ruta:
        return ruta

    # 2. Rutas alternativas obligatorias en distribuciones Ubuntu/Linux
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
    """Concatena todos los modelos HMM individuales en una base de datos única

    y la prensa con hmmpress para permitir el uso de hmmscan.
    """
    db_combinada = os.path.join(carpeta_pfam, "pfam_db.hmm")
    
    archivos_hmm = [
        os.path.join(carpeta_pfam, f) 
        for f in os.listdir(carpeta_pfam) 
        if f.endswith("") and f != "pfam_db.hmm"
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


def generar_logos_hmm(carpeta_pfam, carpeta_salida, ruta_hmmconvert):
    """Extrae los valores de emisión y alturas para la creación de logos

    de secuencias de los perfiles HMM usando hmmconvert de HMMER.
    """
    print("\nGenerando matrices de datos de logotipos para las familias Pfam...")
    archivos_hmm = [
        os.path.join(carpeta_pfam, f) 
        for f in os.listdir(carpeta_pfam) 
        if f.endswith(".hmm") and f != "pfam_db.hmm"
    ]
    
    carpeta_logos = os.path.join(carpeta_salida, "logos_hmm")
    os.makedirs(carpeta_logos, exist_ok=True)
    
    for ruta_hmm in archivos_hmm:
        nombre_familia = os.path.splitext(os.path.basename(ruta_hmm))[0]
        ruta_logo_salida = os.path.join(carpeta_logos, f"{nombre_familia}_logo.txt")
        
        comando = [ruta_hmmconvert, "-w", ruta_hmm]
        
        try:
            with open(ruta_logo_salida, "w", encoding="utf-8") as f_salida:
                subprocess.run(comando, check=True, stdout=f_salida, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print(f"   No se pudo generar el logo para {nombre_familia}")
            
    print(f"¡Matrices de logotipos guardadas con éxito en '{carpeta_logos}'!")


if __name__ == "__main__":
    # Localizar de manera robusta los binarios de HMMER en Ubuntu WSL
    bin_hmmscan = buscar_binario_wsl("hmmscan")
    bin_hmmpress = buscar_binario_wsl("hmmpress")
    bin_hmmconvert = buscar_binario_wsl("hmmconvert")

    if not bin_hmmscan or not bin_hmmpress or not bin_hmmconvert:
        print("\n[ERROR CRÍTICO]: No se encontraron las herramientas de HMMER dentro de WSL Ubuntu.")
        print("Por favor, abre tu terminal de WSL Ubuntu y ejecuta:")
        print("   sudo apt-get update && sudo apt-get install -y hmmer\n")
    else:
        print(f"HMMER detectado en WSL de manera exitosa.")
        print(f"-> Usando: {bin_hmmscan}")
        
        inicializar_entorno()
        try:
            # 1. Preparar e indexar la carpeta de Pfam
            ruta_db_lista = preparar_base_datos_pfam(CARPETA_PFAM, bin_hmmpress)
            
            if ruta_db_lista:
                # 2. Identificar familias correlacionando con la carpeta FASTA
                ejecutar_hmmscan(ruta_db_lista, CARPETA_FASTA, CARPETA_RESULTADOS, bin_hmmscan)
                
                # 3. Extraer los datos del logotipo de secuencia por familia
                generar_logos_hmm(CARPETA_PFAM, CARPETA_RESULTADOS, bin_hmmconvert)
                
                print("\n Pipeline completado con éxito en WSL.")
                print(f"Resultados guardados en la carpeta: '{CARPETA_RESULTADOS}'")
        except Exception as e:
            print(f"\nSe detuvo el proceso debido a un error: {e}")
