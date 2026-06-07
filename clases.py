import os
import re

class ProteinaEscaneada:
    """Clase que representa una proteína y almacena sus hits, tablas crudas y alineamientos completos."""
    def __init__(self, nombre_id):
        self.id_uniprot = nombre_id
        self.hits = []       
        self.reporte_completo = "No se encontró archivo de alineamiento .txt"  
        self.tabla_cruda = "No se encontró archivo de salida tabular .tbl"

    def añadir_hit(self, datos_hit):
        self.hits.append(datos_hit)


class ResultadoHMMER:
    """Clase contenedora principal que agrupa, parsea y exporta todo el pipeline."""
    def __init__(self, carpeta_resultados="resultados_hmmer", carpeta_modelos_origen="modelos_hmm_pfam"):
        self.carpeta = carpeta_resultados
        self.carpeta_modelos = carpeta_modelos_origen
        self.proteinas = {}  
        self._cargar_y_parsear_todo()

    def _limpiar_id_uniprot(self, nombre_archivo) -> str:
        """Extrae de forma robusta el identificador de UniProt."""
        nombre_limpio = nombre_archivo.replace("tabla_", "").replace("resultado_", "")
        match = re.search(r'([A-Z0-9]{6,10})', nombre_limpio)
        if match:
            return match.group(1)
        return os.path.splitext(nombre_limpio)[0]

    def _cargar_y_parsear_todo(self):
        """Carga y mapea de forma cruzada tablas y reportes de alineamientos."""
        if not os.path.exists(self.carpeta):
            print(f"[Error] La carpeta de resultados '{self.carpeta}' no existe.")
            return

        # 1. Parsear tablas (.tbl) de HMMER
        for archivo in os.listdir(self.carpeta):
            if archivo.startswith("tabla_") and archivo.endswith(".tbl"):
                id_prot = self._limpiar_id_uniprot(archivo)
                if id_prot not in self.proteinas:
                    self.proteinas[id_prot] = ProteinaEscaneada(id_prot)
                
                ruta_tbl = os.path.join(self.carpeta, archivo)
                with open(ruta_tbl, "r", encoding="utf-8") as f:
                    self.proteinas[id_prot].tabla_cruda = f.read()
                self._parsear_archivo_tabla(ruta_tbl, self.proteinas[id_prot])

        # 2. Sincronizar reportes (.txt) de alineamientos
        for archivo in os.listdir(self.carpeta):
            if archivo.startswith("resultado_") and archivo.endswith(".txt"):
                id_prot = self._limpiar_id_uniprot(archivo)
                if id_prot not in self.proteinas:
                    self.proteinas[id_prot] = ProteinaEscaneada(id_prot)
                with open(os.path.join(self.carpeta, archivo), "r", encoding="utf-8") as f:
                    self.proteinas[id_prot].reporte_completo = f.read()

    def _parsear_archivo_tabla(self, ruta_archivo, objeto_proteina):
        """Parsea de manera robusta las columnas de un archivo tabular de HMMER."""
        with open(ruta_archivo, "r", encoding="utf-8") as f:
            for linea in f:
                if linea.startswith("#"):
                    continue
                partes = linea.split()
                if len(partes) >= 22:
                    # Limpiar los números de versión agregados por HMMER (.24, .32, etc.)
                    fam_limpia = partes[0].split(".")[0].strip()
                    acc_limpio = partes[1].split(".")[0].strip()
                    
                    hit_info = {
                        "target_name": fam_limpia,      
                        "target_accession": acc_limpio, 
                        "e_value": partes[4],          
                        "score": float(partes[5]),     
                        "ali_from": int(partes[17]),   
                        "ali_to": int(partes[18]),     
                    }
                    objeto_proteina.añadir_hit(hit_info)

    def exportar_resumen_html(self, nombre_salida="resumen_ejecutivo_hmmer.html"):
        """Genera el Dashboard interactivo con el mapeo de secuencias y tablas de HMMER."""
        print(f"-> Ensamblando interfaz web interactiva completa para {len(self.proteinas)} proteínas...")
        
        html_plantilla = """<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <title>Dashboard Analítico HMMER Nativo</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; background-color: #f0f2f5; color: #2d3748; margin: 0; padding: 20px; }}
        .container {{ max-width: 1300px; margin: auto; background: white; padding: 30px; border-radius: 12px; box-shadow: 0 4px 20px rgba(0,0,0,0.08); }}
        h1 {{ color: #1a365d; border-bottom: 3px solid #3182ce; padding-bottom: 10px; margin-top: 0; }}
        
        .collapsible {{ background-color: #2b6cb0; color: white; cursor: pointer; padding: 16px; width: 100%; border: none; text-align: left; outline: none; font-size: 16px; font-weight: bold; margin-top: 10px; border-radius: 6px; display: flex; justify-content: space-between; align-items: center; transition: 0.2s; }}
        .active, .collapsible:hover {{ background-color: #1a365d; }}
        .collapsible:after {{ content: '\\002B'; font-size: 20px; font-weight: bold; float: right; margin-left: 5px; }}
        .active:after {{ content: "\\2212"; }}
        .content {{ padding: 0 18px; display: none; overflow: hidden; background-color: #fff; border: 1px solid #e2e8f0; border-top: none; border-bottom-left-radius: 6px; border-bottom-right-radius: 6px; }}
        
        .tab {{ overflow: hidden; border-bottom: 2px solid #e2e8f0; margin-top: 15px; display: flex; gap: 5px; }}
        .tab button {{ background-color: #edf2f7; color: #4a5568; border: none; outline: none; cursor: pointer; padding: 10px 20px; font-size: 14px; font-weight: 600; border-top-left-radius: 6px; border-top-right-radius: 6px; transition: 0.2s; }}
        .tab button:hover {{ background-color: #e2e8f0; }}
        .tab button.active-sub {{ background-color: #3182ce; color: white; }}
        .tabcontent {{ display: none; padding: 15px 0; }}
        
        pre {{ background-color: #1a202c; color: #edf2f7; padding: 15px; border-radius: 6px; overflow-x: auto; font-family: 'Courier New', Courier, monospace; font-size: 13px; max-height: 400px; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
        th {{ background-color: #4a5568; color: white; padding: 10px; text-align: left; font-size: 14px; }}
        td {{ padding: 10px; border-bottom: 1px solid #e2e8f0; font-size: 13px; }}
        
        .badge {{ background: #e2e8f0; color: #2d3748; padding: 3px 8px; border-radius: 4px; font-family: monospace; font-weight: bold; }}
        .badge-ev {{ background: #fed7d7; color: #9b2c2c; padding: 3px 8px; border-radius: 12px; font-weight: bold; }}
        footer {{ text-align: center; margin-top: 40px; font-size: 12px; color: #718096; border-top: 1px solid #e2e8f0; padding-top: 15px; }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 Dashboard de Identificación de Dominios (Mapeo de Hits)</h1>
    <p>Haz clic sobre cualquiera de las secuencias de proteínas mapeadas a continuación para inspeccionar sus alineamientos detallados y tablas de asignación de familias.</p>
"""
        cuerpo_html = ""
        for index, (id_prot, proteina) in enumerate(sorted(self.proteinas.items())):
            total_hits = len(proteina.hits)
            
            cuerpo_html += f"""
    <button class="collapsible">Proteína / Secuencia: {id_prot} <span style="font-size:13px; background:#fff; color:#2b6cb0; padding:2px 8px; border-radius:10px;">Hits: {total_hits}</span></button>
    <div class="content">
        <div class="tab">
            <button class="tablinks active-sub" onclick="abrirSubPestaña(event, 'hits_{index}')">Tabla de Hits (.tbl)</button>
            <button class="tablinks" onclick="abrirSubPestaña(event, 'ali_{index}')">Alineamiento Completo (.txt)</button>
        </div>
        
        <!-- PESTAÑA 1: TABLA ESTRUCTURADA Y TEXTO CRUDO TBL -->
        <div id="hits_{index}" class="tabcontent" style="display:block;">
            <h3>Resultados tabulares filtrados de HMMER</h3>
            <table>
                <thead>
                    <tr>
                        <th>Familia (NAME)</th>
                        <th>Acceso (ACC)</th>
                        <th>E-Value Global</th>
                        <th>Score (bits)</th>
                        <th>Coordenadas Dominio</th>
                    </tr>
                </thead>
                <tbody>"""
                
            if not proteina.hits:
                cuerpo_html += "<tr><td colspan='5' style='color:#a0aec0; font-style:italic;'>No se detectó ningún dominio homólogo con significancia estadística.</td></tr>"
            else:
                for hit in proteina.hits:
                    cuerpo_html += f"""
                    <tr>
                        <td><span class="badge">{hit['target_name']}</span></td>
                        <td><span class="badge">{hit['target_accession']}</span></td>
                        <td><span class="badge-ev">{hit['e_value']}</span></td>
                        <td>{hit['score']} b</td>
                        <td>Residuos: {hit['ali_from']} - {hit['ali_to']}</td>
                    </tr>"""
            
            cuerpo_html += f"""
                </tbody>
            </table>
            <h4>Texto crudo original extraído del archivo tabular (.tbl):</h4>
            <pre>{proteina.tabla_cruda}</pre>
        </div>
        
        <!-- PESTAÑA 2: ALINEAMIENTO COMPLETO TXT -->
        <div id="ali_{index}" class="tabcontent">
            <h3>Alineamiento de Homología de Markov por Posición</h3>
            <pre>{proteina.reporte_completo}</pre>
        </div>
    </div>"""

        script_js = """
    <footer>Reporte de Alineamientos Integrado • Ejecutado bajo el entorno WSL Ubuntu</footer>
</div>

<script>
    var coll = document.getElementsByClassName("collapsible");
    for (var i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display === "block") {
                content.style.display = "none";
            } else {
                content.style.display = "block";
            }
        });
    }

    // Lógica para alternar las sub-pestañas internas de datos crudos o tablas
    function abrirSubPestaña(evt, nombrePestaña) {
        var tabcontent, tablinks;
        var contenedorPadre = evt.currentTarget.parentElement.parentElement;
        
        tabcontent = contenedorPadre.getElementsByClassName("tabcontent");
        for (var i = 0; i < tabcontent.length; i++) {
            tabcontent[i].style.display = "none";
        }
        
        tablinks = contenedorPadre.getElementsByClassName("tablinks");
        for (var i = 0; i < tablinks.length; i++) {
            tablinks[i].classList.remove("active-sub");
        }
        
        document.getElementById(nombrePestaña).style.display = "block";
        evt.currentTarget.classList.add("active-sub");
    }
</script>
</body>
</html>
"""
        # Unificación de los bloques para compilar el documento HTML final
        html_final = html_plantilla + cuerpo_html + script_js
        
        with open(nombre_salida, "w", encoding="utf-8") as f_html:
            f_html.write(html_final)
        print(f"¡Éxito! Reporte estructurado guardado con éxito para {len(self.proteinas)} proteínas.")


# PUNTO DE ENTRADA CORREGIDO CON SINTAXIS NATIVA DE PYTHON (DOBLE GUION BAJO)
if __name__ == "__main__":
    datos = ResultadoHMMER("resultados_hmmer")
    datos.exportar_resumen_html()
