import os
import re


class FamiliaLogo:
    """Clase que almacena la matriz de información de un logotipo generado por hmmlogo."""

    def __init__(self, ruta_logo):
        self.id_familia = os.path.basename(ruta_logo).replace("_logo.txt", "")
        self.posiciones = {}  # Diccionario: {num_posicion: {aminoacido: bits}}
        self._parsear_archivo_logo(ruta_logo)

    def _parsear_archivo_logo(self, ruta):
        if not os.path.exists(ruta):
            return
        
        with open(ruta, "r", encoding="utf-8") as f:
            lineas = f.readlines()
            
        # Buscar la cabecera de aminoácidos
        aminoacidos = []
        for linea in lineas:
            if linea.startswith("pos") or "m->i" in linea:
                continue
            # Línea de cabecera de aminoácidos (ej: "     A      C      D...")
            if len(aminoacidos) == 0 and any(aa in linea for aa in ["A", "C", "D", "E"]):
                aminoacidos = linea.split()
                continue
            
            # Líneas de datos (ej: "   1   0.12   0.00   1.45...")
            partes = linea.split()
            if partes and partes[0].isdigit() and len(partes) > len(aminoacidos):
                pos = int(partes[0])
                valores = [float(x) for x in partes[1:len(aminoacidos)+1]]
                self.posiciones[pos] = dict(zip(aminoacidos, valores))

    def obtener_residuos_principales(self, posicion, top_n=3):
        """Devuelve los N aminoácidos más conservados en una posición específica."""
        if posicion not in self.posiciones:
            return []
        residuos_ordenados = sorted(
            self.posiciones[posicion].items(), key=lambda item: item[1], reverse=True
        )
        return [res for res in residuos_ordenados if res[1] > 0][:top_n]


class ProteinaEscaneada:
    """Clase que representa una proteína y almacena sus hits/familias identificadas."""

    def __init__(self, nombre_id):
        self.id_uniprot = nombre_id
        self.hits = []       # Lista de diccionarios con información de dominios detectados
        self.reporte_completo = ""  # Texto crudo del archivo de alineamiento .txt

    def añadir_hit(self, datos_hit):
        """Añade un hit estructurado proveniente del archivo .tbl."""
        self.hits.append(datos_hit)

    def obtener_familias_validas(self, e_value_corte=1e-5):
        """Filtra y devuelve los nombres de las familias que superan un umbral de E-value."""
        return [
            hit["target_name"] 
            for hit in self.hits 
            if float(hit["e_value"]) <= e_value_corte
        ]


class ResultadoHMMER:
    """Clase contenedora principal (Master Class) que agrupa y parsea todo

    el directorio 'resultados_hmmer'.
    """

    def __init__(self, carpeta_resultados="resultados_hmmer"):
        self.carpeta = carpeta_resultados
        self.proteinas = {}  # {id_uniprot: Objeto ProteinaEscaneada}
        self.logos = {}      # {id_familia: Objeto FamiliaLogo}
        
        # Ejecutar el parseo completo automáticamente al instanciar la clase
        self._cargar_y_parsear_todo()

    def _cargar_y_parsear_todo(self):
        if not os.path.exists(self.carpeta):
            print(f"[Aviso] La carpeta '{self.carpeta}' no existe todavía.")
            return

        # 1. Parsear archivos tabulares (.tbl) de hmmscan
        for archivo in os.listdir(self.carpeta):
            if archivo.startswith("tabla_") and archivo.endswith(".tbl"):
                # Extrae el ID de la proteína eliminando 'tabla_' y '.tbl'
                id_prot = archivo.replace("tabla_", "").split(".")[0]
                
                if id_prot not in self.proteinas:
                    self.proteinas[id_prot] = ProteinaEscaneada(id_prot)
                
                self._parsear_archivo_tabla(os.path.join(self.carpeta, archivo), self.proteinas[id_prot])

            # 2. Guardar reportes de texto planos (.txt) de alineamientos
            elif archivo.startswith("resultado_") and archivo.endswith(".txt"):
                id_prot = archivo.replace("resultado_", "").split(".")[0]
                
                if id_prot not in self.proteinas:
                    self.proteinas[id_prot] = ProteinaEscaneada(id_prot)
                    
                with open(os.path.join(self.carpeta, archivo), "r", encoding="utf-8") as f:
                    self.proteinas[id_prot].reporte_completo = f.read()

        # 3. Parsear subcarpeta de logos_hmm si existe
        ruta_logos = os.path.join(self.carpeta, "logos_hmm")
        if os.path.exists(ruta_logos):
            for archivo in os.listdir(ruta_logos):
                if archivo.endswith("_logo.txt"):
                    objeto_logo = FamiliaLogo(os.path.join(ruta_logos, archivo))
                    self.logos[objeto_logo.id_familia] = objeto_logo

    def _parsear_archivo_tabla(self, ruta_archivo, objeto_proteina):
        """Parsea de manera robusta las columnas de un archivo --domtblout de HMMER."""
        with open(ruta_archivo, "r", encoding="utf-8") as f:
            for linea in f:
                if linea.startswith("#"):  # Saltar comentarios de HMMER
                    continue
                partes = linea.split()
                if len(partes) >= 22:
                    # Estructuración de columnas esenciales de HMMER según su documentación oficial
                    hit_info = {
                        "target_name": partes[0],      # Familia de Pfam (Ej: PF00069)
                        "target_accession": partes[1], # Acceso Pfam alterno
                        "e_value": partes[6],          # E-value global de la secuencia
                        "score": float(partes[7]),     # Puntuación en bits de HMMER
                        "ali_from": int(partes[17]),   # Coordenada inicio de alineamiento en proteína
                        "ali_to": int(partes[18]),     # Coordenada fin de alineamiento en proteína
                    }
                    objeto_proteina.añadir_hit(hit_info)

