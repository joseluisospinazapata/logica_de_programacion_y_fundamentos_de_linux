import os

class FamiliaLogo:
    """Clase que almacena la matriz de información de un logotipo generado por hmmlogo."""
    def __init__(self, ruta_logo):
        self.id_familia = os.path.basename(ruta_logo).replace("_logo.txt", "")
        self.posiciones = {}  
        self._parsear_archivo_logo(ruta_logo)

    def _parsear_archivo_logo(self, ruta):
        if not os.path.exists(ruta):
            return
        with open(ruta, "r", encoding="utf-8") as f:
            lineas = f.readlines()
            
        aminoacidos = []
        for linea in lineas:
            # Saltar líneas informativas iniciales o de transición de gaps
            if linea.startswith("pos") or "m->i" in linea or linea.strip().startswith("#"):
                continue
            
            partes = linea.split()
            if not partes:
                continue
                
            # Identificar la cabecera de letras de los aminoácidos (ej: "A C D E F G...")
            if len(aminoacidos) == 0 and any(aa in partes for aa in ["A", "C", "D", "E"]):
                aminoacidos = partes
                continue
            
            # CORRECCIÓN: Comprobar si el primer elemento de la lista (el índice de posición) es un número
            if partes[0].isdigit() and len(partes) > len(aminoacidos):
                pos = int(partes[0])
                # Mapear los valores numéricos de bits a cada aminoácido correspondiente
                valores = [float(x) for x in partes[1:len(aminoacidos)+1]]
                self.posiciones[pos] = dict(zip(aminoacidos, valores))

    def obtener_residuos_principales(self, posicion, top_n=3):
        """Devuelve los N aminoácidos más conservados en una posición específica."""
        if posicion not in self.posiciones:
            return []
        # Ordenar los residuos de mayor a menor puntuación en bits
        residuos_ordenados = sorted(
            self.posiciones[posicion].items(), key=lambda item: item[1], reverse=True
        )
        return [f"{res} ({val:.2f} b)" for res, val in residuos_ordenados if val > 0][:top_n]


class ProteinaEscaneada:
    """Clase que representa una proteína y almacena sus hits/familias identificadas."""
    def __init__(self, nombre_id):
        self.id_uniprot = nombre_id
        self.hits = []       
        self.reporte_completo = ""  

    def añadir_hit(self, datos_hit):
        self.hits.append(datos_hit)

    def obtener_familias_validas(self, e_value_corte=1e-5):
        return [hit["target_name"] for hit in self.hits if float(hit["e_value"]) <= e_value_corte]


class ResultadoHMMER:
    """Clase contenedora principal que agrupa, parsea y exporta los datos."""
    def __init__(self, carpeta_resultados="resultados_hmmer"):
        self.carpeta = carpeta_resultados
        self.proteinas = {}  
        self.logos = {}      
        self._cargar_y_parsear_todo()

    def _cargar_y_parsear_todo(self):
        if not os.path.exists(self.carpeta):
            return

        for archivo in os.listdir(self.carpeta):
            if archivo.startswith("tabla_") and archivo.endswith(".tbl"):
                id_prot = archivo.replace("tabla_", "").replace(".tbl", "")
                if id_prot not in self.proteinas:
                    self.proteinas[id_prot] = ProteinaEscaneada(id_prot)
                self._parsear_archivo_tabla(os.path.join(self.carpeta, archivo), self.proteinas[id_prot])

            elif archivo.startswith("resultado_") and archivo.endswith(".txt"):
                id_prot = archivo.replace("resultado_", "").replace(".txt", "")
                if id_prot not in self.proteinas:
                    self.proteinas[id_prot] = ProteinaEscaneada(id_prot)
                with open(os.path.join(self.carpeta, archivo), "r", encoding="utf-8") as f:
                    self.proteinas[id_prot].reporte_completo = f.read()

        ruta_logos = os.path.join(self.carpeta, "logos_hmm")
        if os.path.exists(ruta_logos):
            for archivo in os.listdir(ruta_logos):
                if archivo.endswith("_logo.txt"):
                    objeto_logo = FamiliaLogo(os.path.join(ruta_logos, archivo))
                    self.logos[objeto_logo.id_familia] = objeto_logo

    def _parsear_archivo_tabla(self, ruta_archivo, objeto_proteina):
        with open(ruta_archivo, "r", encoding="utf-8") as f:
            for linea in f:
                if linea.startswith("#"):
                    continue
                partes = linea.split()
                if len(partes) >= 22:
                    hit_info = {
                        "target_name": partes[0],      
                        "target_accession": partes[1], 
                        "e_value": partes[6],          
                        "score": float(partes[7]),     
                        "ali_from": int(partes[17]),   
                        "ali_to": int(partes[18]),     
                    }
                    objeto_proteina.añadir_hit(hit_info)

    def exportar_resumen_html(self, nombre_salida="resumen_ejecutivo_hmmer.html"):
        """Genera un reporte ejecutivo en HTML con diseño responsivo y moderno."""
        print(f"-> Generando reporte ejecutivo en HTML: '{nombre_salida}'...")
        
        html_plantilla = """<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Reporte Ejecutivo - Pipeline HMMER & UniProt</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; background-color: #f4f6f9; color: #333; margin: 0; padding: 30px; }}
        .container {{ max-width: 1200px; margin: auto; background: white; padding: 30px; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); }}
        h1 {{ color: #1a365d; border-bottom: 3px solid #2b6cb0; padding-bottom: 10px; margin-top: 0; }}
        .stats {{ display: flex; gap: 20px; margin-bottom: 25px; }}
        .card {{ background: #ebf8ff; border: 1px solid #bee3f8; padding: 15px 25px; border-radius: 8px; flex: 1; text-align: center; }}
        .card h3 {{ margin: 0; color: #2c5282; font-size: 14px; text-transform: uppercase; }}
        .card p {{ margin: 5px 0 0 0; font-size: 28px; font-weight: bold; color: #2b6cb0; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 15px; background: white; }}
        th {{ background-color: #2b6cb0; color: white; text-align: left; padding: 12px; font-weight: 600; }}
        td {{ padding: 12px; border-bottom: 1px solid #e2e8f0; font-size: 14px; }}
        tr:hover {{ background-color: #f7fafc; }}
        .badge-id {{ background: #edf2f7; padding: 4px 8px; border-radius: 4px; font-family: monospace; font-weight: bold; color: #4a5568; }}
        .badge-fam {{ background: #e2e8f0; padding: 4px 8px; border-radius: 4px; font-family: monospace; color: #2d3748; font-weight: 500; }}
        .badge-ev {{ background: #c6f6d5; color: #22543d; padding: 3px 8px; border-radius: 20px; font-size: 12px; font-weight: bold; }}
        .logo-box {{ background: #fffaf0; border: 1px solid #feebc8; padding: 6px; border-radius: 6px; font-size: 12px; font-family: monospace; color: #dd6b20; }}
        footer {{ margin-top: 40px; text-align: center; font-size: 12px; color: #a0aec0; padding-top: 20px; border-top: 1px solid #e2e8f0; }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 Resumen Ejecutivo: Mapeo de Familias Proteicas</h1>
    <div class="stats">
        <div class="card"><h3>Proteínas Analizadas</h3><p>{total_prot}</p></div>
        <div class="card"><h3>Familias Pfam Identificadas</h3><p>{total_fam}</p></div>
    </div>
    <table>
        <thead>
            <tr>
                <th>ID UniProt</th>
                <th>Familia Identificada (Hit)</th>
                <th>E-Value</th>
                <th>Score (bits)</th>
                <th>Rango de Alineamiento</th>
                <th>Residuos Conservados del Logo (Muestra)</th>
            </tr>
        </thead>
        <tbody>
"""
        filas = ""
        familias_unicas = set()
        
        for id_prot, proteina in sorted(self.proteinas.items()):
            if not proteina.hits:
                filas += f"<tr><td><span class='badge-id'>{id_prot}</span></td><td colspan='5' style='color:#a0aec0; font-style:italic;'>Sin dominios identificados</td></tr>"
                continue
                
            for hit in proteina.hits:
                fam = hit["target_name"]
                familias_unicas.add(fam)
                
                logo_str = "Matriz N/A"
                if fam in self.logos:
                    logo_obj = self.logos[fam]
                    # Extraer el residuo dominante principal de las primeras 3 posiciones del dominio
                    pos1 = logo_obj.obtener_residuos_principales(1, top_n=1)
                    pos2 = logo_obj.obtener_residuos_principales(2, top_n=1)
                    pos3 = logo_obj.obtener_residuos_principales(3, top_n=1)
                    
                    p1 = pos1[0] if pos1 else "-"
                    p2 = pos2[0] if pos2 else "-"
                    p3 = pos3[0] if pos3 else "-"
                    logo_str = f"Pos1: {p1} | Pos2: {p2} | Pos3: {p3}"

                filas += f"""
            <tr>
                <td><span class='badge-id'>{id_prot}</span></td>
                <td><span class='badge-fam'>{fam}</span></td>
                <td><span class='badge-ev'>{hit['e_value']}</span></td>
                <td>{hit['score']}</td>
                <td>AA: {hit['ali_from']} - {hit['ali_to']}</td>
                <td><div class='logo-box'>{logo_str}</div></td>
            </tr>"""

        html_final = html_plantilla.format(
            total_prot=len(self.proteinas),
            total_fam=len(familias_unicas)
        ) + filas + """
        </tbody>
    </table>
    <footer>Reporte generado automáticamente por el Pipeline Científico Integrado</footer>
</div>
</body>
</html>
"""
        with open(nombre_salida, "w", encoding="utf-8") as f_html:
            f_html.write(html_final)

print(f"¡Éxito! Reporte web guardado en '{os.path.abspath('resumen_ejecutivo_hmmer.html')}'.")

if __name__ == "main":
    datos = ResultadoHMMER("resultados_hmmer")
    datos.exportar_resumen_html()