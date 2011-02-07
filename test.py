from GOGraph import *

location = "./go_daily-termdb.obo-xml"
#location = "./testgraph.obo-xml"

#possible namespaces are: "biological_process", "molecular_function", "cellular_component"
GOGraphtest = GOGraph("biological_process", location)
