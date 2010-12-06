from GOGraph import *

location = "./go_daily-termdb.obo-xml"

#possible namespaces are: "biological_process", "molecular_function", "cellular_component"
GOGraph = GOGraph("biological_process", location)
