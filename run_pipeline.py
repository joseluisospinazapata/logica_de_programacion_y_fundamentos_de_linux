import os
import sys
import subprocess

# Configuración de los módulos en orden secuencial de ejecución
MODULOS_PIPELINE = [
    "descarga_pfam.py",            # Fase 1: Descarga modelos HMM de Pfam
    "descarga_UniProt.py",         # Fase 2: Descarga archivos FASTA y gestiona ETags
    "alineamientos_logos_hmmer.py" # Fase 3: Ejecuta hmmpress, hmmscan, hmmlogo y limpia binarios
]


def ejecutar_script_python(nombre_script):
    """Ejecuta un script de Python de forma controlada y transmite su salida en vivo."""
    print("\n" + "="*70)
    print(f" EJECUTANDO MÓDULO: {nombre_script}")
    print("="*70)
    
    if not os.path.exists(nombre_script):
        print(f"[ERROR CRÍTICO]: No se encuentra el archivo '{nombre_script}' en el directorio.")
        return False

    # Ejecuta forzando el uso del intérprete de Python activo en WSL
    resultado = subprocess.run([sys.executable, nombre_script])
    
    if resultado.returncode == 0:
        print(f"-> Módulo {nombre_script} finalizado correctamente.")
        return True
    else:
        print(f"[ERROR]: El script '{nombre_script}' falló con código de salida: {resultado.returncode}")
        return False


def lanzar_pipeline_completo():
    print("======================================================================")
    print("        SISTEMA DE CONTROL DE PIPELINE BIOINFORMÁTICO INTEGRADO       ")
    print("======================================================================")

    # 1. Ejecutar de forma encadenada los 3 módulos analíticos previos
    for script in MODULOS_PIPELINE:
        exito = ejecutar_script_python(script)
        if not exito:
            print(f"\n[PIPELINE ABORTADO] Se detuvo el flujo debido a un fallo en: {script}")
            sys.exit(1)

    print("\n" + "="*70)
    print(" FASE FINAL: INSTANCIACIÓN DE CLASES Y REPORTE EJECUTIVO")
    print("="*70)

    # 2. Invocar e instanciar dinámicamente el generador HTML de la clase principal
    try:
        # Importación dinámica del módulo de clases que acabamos de guardar
        from clases import ResultadoHMMER
        
        print("-> Cargando e indexando la carpeta de resultados en la memoria orientada a objetos...")
        instancia_master = ResultadoHMMER("resultados_hmmer")
        
        # Exportación del reporte ejecutivo HTML final
        instancia_master.exportar_resumen_html("resumen_ejecutivo_hmmer.html")
        
        print("\n======================================================================")
        print(" 🎉 [PROCESO COMPLETADO] ¡Fases analíticas y reporte web finalizados! ")
        print("======================================================================")
        
    except ImportError:
        print("[ERROR FASE FINAL]: No se pudo importar 'clases.py'. Asegúrate de que el archivo exista.")
    except Exception as e:
        print(f"[ERROR FASE FINAL]: Ocurrió un fallo al compilar los objetos y el HTML: {e}")


if __name__ == "__main__":
    lanzar_pipeline_completo()
