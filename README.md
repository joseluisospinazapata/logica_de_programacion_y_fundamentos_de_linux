# Materia: logica de programacion y fundamentos de linux
Repositorio sobre el curso lógica de programación y fundamentos en linux del post grado en bioinformática institución universitaria ITM primera cohorte 2026.

# Trabajo final del curso

# Pipeline bioinformático para identificación de dominios de familias de proteínas en secuencias de aminoácidos

El siguiente proceso describe el desarrollo de un problema formulado para el curso en donde por etapas en lenguaje de programación Python o R optativos, se ejecutara a travez de un pipeline con resultados finales en formato HTML, la identificación de dominios con HMMER de una serie de familias de proteínas Pfam en una serie de secuencias de aminoácidos identificadores UniProt.

# Familias Pfam: 
Protein_kinase, Pkinase_Tyr, Ras, SH2, SH3_1, zf-C2H2, Homeobox, HTH_1, bZIP_1, Myb_DNA-binding, RRM_1, DEAD, KH_1, dsrm, ABC_tran, MFS_1, Ion_trans, HlyD, Aminotran_1_2, Aldedh, TIM, NAD_binding_1, GST_C_family, WD40, Ank, TPR_1, LRR_1, HSP70, HSP20, DnaJ, Response_reg, HisKA, Peptidase_M16, Sigma70_r2, Immunoglobulin, EGF, Cadherin, Fibronectin.

# Identificadores UniProt:
P00519, P42684, P12931, P06241, P07947, Q06187, P43403, P43405, P62993, P01112, P01116, P01111, P04049, P31749, P28482, P27361, P00533, P21802, P16234, P12956, P29353, P42681, P35222, P62937, P29317, P08047, P15056, P40763, P42224, P15924, P10242, P19838, P11473, P61244, Q9Y2T1, P08107, P0A6Y8, P0A6W5, P0A9Q7, P0A799, P0A7Y4, P0A8V2, P39451, P11142, P13569, P22681, P98160, P12814, Q92793, Q13485.

# Diagrama de flujo:
<img width="721" height="463" alt="image" src="https://github.com/user-attachments/assets/52d1af1d-faf2-4d2b-998b-029c8ef26c86" />

# Requerimientos del sistema:

# 1. Entorno Base (WSL, VS Code y Python)

# 1.1 Windows Subsystem for Linux (WSL) y Ubuntu
Como la herramienta central (HMMER) está hecha de forma nativa para Linux, prepararemos un entorno Ubuntu dentro de tu Windows.
Instalación: Abre la terminal PowerShell de Windows como administrador y escribe el comando: "wsl --install -d Ubuntu".
Reinicia la computadora si el sistema lo pide. Al abrir Ubuntu por primera vez, te solicitará crear un usuario y una contraseña. Luego abre la aplicación "Ubuntu" en el menú de inicio para interactuar con la terminal Linux mediante comandos de consola.

# 1.2 Visual Studio Code (VS Code)
Instalación: Descarga e instala el asistente de Windows desde la Página Oficial de VS Code. Una vez dentro de VS Code, ve a la sección de Extensiones (icono de cuadrados a la izquierda) e instala la extensión llamada WSL. Luego en la terminal de Ubuntu escribe "code .". Esto abrirá VS Code de forma automática en Windows, pero ejecutando los archivos y las herramientas directamente dentro del sistema Linux de forma integrada.

# 1.3 Python
Instalación: Ubuntu ya incluye Python por defecto. Solo se necesita instalar su gestor de paquetes (pip) en la terminal de Ubuntu corriendo el comando: "sudo apt update && sudo apt install python3-pip python3-venv -y", Con esta herramienta se crearán los modulos con extensión .py dentro de VS Code para gestionar y enlazar los archivos del pipeline.

# 2. Bases de Datos Biológicas (Entradas)
No se requiere instalar programas para las bases de datos de entrada; Se descargaran por medio de modulos especificos de Python a travez de sus correspondientes APIs, las secuencias y matrices en carpetas locales del proyecto en los formatos estándar de bioinformática.

# 3. HMMER (Versión ejecutable en local)
Instalación: En la terminal de Ubuntu de WSL, puedes instalar HMMER directamente con el gestor de paquetes de Linux ejecutando el comando:
"sudo apt install hmmer -y"

# 4. Formato compatible con HTML5
Un modulo especifico de Python traducirá los resultados planos de HMMER a un archivo web visual (HTML). Para abrirlos solo necesitas cualquier navegador web moderno (como Google Chrome, Firefox o Microsoft Edge).

# Etapas de la programacion del pipeline

# Etapa 1 descarga.

Inicialmente se crearon dos modulos separados para esta etapa con el proposito de diferenciar las bases de datos oficiales de las cuales se sustraeran los datos especificos del problema.

# Modulo 1 descarga_pfam.py

Este código sirve para descargar de forma automática los modelos de las familias de proteínas desde el sitio web oficial de Pfam con los id. accesion de cada familia modelo. El script busca los modelos de una lista de 38 id. accesion de familias de proteínas Pfam, guarda la información en la particion de linux de computadora y la deja lista para usar.

# Pasos del script:
# 1. Configuración de seguridad
Conexión inteligente: Prepara la descarga para que no falle fácilmente. Si el servidor de internet da un error temporal, el script lo vuelve a intentar hasta 3 veces automáticamente. 
Pausa de cortesía: Espera medio segundo entre cada descarga para no saturar o sobrecargar el servidor web.

#  2. Control de archivos existentes
Evita repetir trabajo: Antes de descargar una familia de proteínas, revisa si el archivo ya existe en la carpeta. Si ya está ahí, se lo salta para ahorrar tiempo.

# 3. Descarga y revisión del formato
Descarga: Conecta con la base de datos mediante una dirección web dinámica para cada familia de proteína.
Verificación: Revisa los primeros datos internos (bytes mágicos) para saber si viene comprimido (formato .gz) o si es texto normal.
Procesamiento: Si está comprimido, lo descomprime. Si es texto normal, solo le cambia el nombre. Al final, todos los archivos quedan guardados con el formato correcto que se necesita (.hmm).

# 4. Registro de actividades (Logs)
Bitácora: Todo lo que pasa se anota en la consola y en un archivo de texto llamado "registro_descarga.log", Escribe la hora exacta de cada descarga, si tuvo éxito, si falló por culpa de internet o si el archivo ya estaba repetido.

# Modulo 2 descarga_UniProt.py

Este código sirve para descargar las secuencias de proteínas en formato FASTA desde la base de datos de UniProt. El script revisa una lista de 50 identificadoes de proteínas UniProt, las guarda en la particion de linux de la computadora y utiliza un sistema inteligente para no repetir descargas a menos que la proteína haya cambiado en la base de datos.

# Pasos del script:
# 1. Control inteligente de versiones (ETags)
No descarga por descargar: Antes de bajar un archivo, el código hace una pregunta rápida al servidor (HEAD) para revisar la "etiqueta de versión" (llamada ETag) de la proteína. 
Compara y decide: Guarda estas etiquetas en un archivo llamado "control_versiones_uniprot.json". Si la proteína no ha cambiado en la base de datos, se salta la descarga para ahorrar tiempo.
Detecta cambios: Si la proteína fue modificada o actualizada en UniProt, el código lo detecta, lanza una alerta y vuelve a descargar la nueva versión.

# 2. Seguridad en las descargas
Reintentos automáticos: Si la conexión a internet falla por un segundo, el script no se rompe. Vuelve a intentar la descarga hasta 3 veces automáticamente.
Pausa de cortesia: Espera 0.3 segundos entre cada proteína para no saturar el servidor de UniProt.

# 3. Verificación de archivos reales
Filtra errores: Cuando descarga el archivo, el script verifica que no esté vacío y que contenga el símbolo ">", que es el inicio obligatorio de cualquier archivo de proteínas tipo FASTA. Si no lo tiene, lo rechaza. 

# 4. Reporte de fallas
Lista negra: Si una proteína da error al descargar, no existe o el archivo está dañado, el código lo anota todo en un archivo de texto llamado "errores_uniprot.log" junto con la hora exacta y el motivo del fallo.

# Etapa 2 procesamiento

# Modulo 3 alineamientos_hmmer.py
Este código sirve para analizar secuencias de proteínas utilizando la herramienta HMMER. El script toma las familias de proteínas (Pfam) y las secuencias (FASTA) descargadas, las compara entre sí y guarda los resultados del análisis en una carpeta nueva.

# Pasos del script:
# 1. Validación de herramientas
Busca programas: El código revisa si tu computadora tiene instalados dos programas esenciales llamados "hmmscan" y "hmmpress" (usando el sistema WSL de Ubuntu).
Filtra errores: Si no encuentra estos programas, el código se detiene y te da las instrucciones exactas en pantalla para instalarlos.

# 2. Creación de una base de datos unificada
Une los archivos: Toma todos los archivos de familias de proteínas (.hmm) individuales descargados en la carpeta "modelos_hmm_pfam" y crea una base de datos temporal llamada "pfam_db.hmm".
Optimiza la lectura: Utiliza la herramienta hmmpress para indexar y preparar esa base de datos unificada. Esto hace que las búsquedas futuras sean muchísimo más rápidas.

# 3. Escaneo y análisis de proteínas
Compara datos: El código toma cada archivo de secuencia de proteína (.fasta) y lo compara contra la base de datos unificada de familias de proteínas.
Genera reportes: Para cada secuencia analizada, crea dos archivos de texto con los resultados detallados dentro de la carpeta "resultados_hmmer" (resultado_...txt y tabla_...tbl).

# 4. Limpieza automática
Borra la basura: Una vez que termina todo el análisis, el script borra de forma segura la base de datos unificada y los archivos temporales binarios que creó al inicio para no ocupar espacio innecesario en tu disco duro.

# Etapa 3 clasificacion

# Modulo 4 clases.py
Este código sirve para organizar los resultados de los análisis de HMMER y crear una página web interactiva (Dashboard) donde revisarlos fácilmente. El script lee los archivos de texto y tablas que generó el programa, extrae los datos importantes y los acomoda en una interfaz visual cómoda.

# Pasos del script:
# 1. Organización de los datos (Clases)
Contenedores de información: Utiliza dos clases (ProteinaEscaneada y ResultadoHMMER). Su trabajo es crear carpetas virtuales en la particion de linux de la computadora para agrupar cada proteína con sus respectivos archivos de texto (.txt) y tablas de datos (.tbl).
Limpieza de nombres: Extrae y limpia de forma automática el identificador único de la proteína (ID de UniProt) a partir del nombre de los archivos, quitando extensiones y texto extra.

# 2. Lectura y extracción (Parsing)
Escaneo de carpetas: Busca dentro de la carpeta de resultados y abre cada archivo guardado.
Filtro de columnas: Lee las tablas de HMMER línea por línea. Ignora los comentarios y extrae únicamente los datos valiosos, como el nombre de la familia de la proteína, el valor de confianza estadística (E-value), la puntuación (Score) y las posiciones exactas de inicio y fin del alineamiento.

# 3. Creación del Dashboard Web
Diseña una página web: Aglomera todos los resultados extraídos y genera de forma automática un archivo interactivo llamado "resumen_ejecutivo_hmmer.html".
Interfaz amigable: La página web creada incluye: Botones desplegables para cada proteína, pestañas para cambiar entre la tabla de resultados y el archivo original con el alineamiento completo en texto plano y por ultimo un diseño visual limpio, tablas y etiquetas para facilitar la lectura.

# Etapa 4 ejecucion

# Modulo 5 pipeline.py
Este modulo es el controlador principal de todo el proyecto bioinformático. Su función es unir los cuatro modulos que se revisaron antes y ejecutarlos uno detrás del otro de forma automática y ordenada.

# Pasos del script:
# 1. Lista de tareas en orden (Pipeline)
Define las fases: El código invoca los tres primeros modulos en el orden que deben correr: Descarga de familias de proteínas (descarga_pfam.py), descarga de secuencias de proteínas (descarga_UniProt.py) y análisis con HMMER (alineamientos_hmmer.py).
Sistema de freno de mano: Si uno de los scripts falla o da un error, el código detiene todo el proceso de inmediato y te avisa cuál fue el falló, esto evita que el análisis continúe con datos incompletos o erróneos.

# 2. Generación del reporte final
Crea la página web: Si los tres modulos anteriores terminan con éxito, el código llama automáticamente al cuarto modulo (clases.py).
Arma el Dashboard: Lee todos los resultados guardados en la particion de linux y activa la función para crear el archivo interactivo resumen_ejecutivo_hmmer.html.

Puedes ver el [Dashboard Interactivo de resumen ejecutivo hmmer aquí](https://joseluisospinazapata.github.io/logica_de_programacion_y_fundamentos_de_linux/resumen_ejecutivo_hmmer.html)











