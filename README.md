# tp-interactomics


## Contexte biologique

Nous allons reproduire certaines analyses rapportées dans une étude du réseau d'interaction des protéines du [virus d'Epstein-Bar](https://en.wikipedia.org/wiki/Epstein%E2%80%93Barr_virus) avec certaines protéines de son hôte principal l'homme ([Calderwood et al.](https://www.pnas.org/content/104/18/7606)). Le cible cellulaire principale du virus EBV est lymphocyte B humain.  

## Mise en place

Seul Jupyter avec la librarire [networkx](https://networkx.org/) est requis.
Les libraries Pandas et request peuvent également être utilisées.
Il est entendu qu'a l'exception des encarts du présent `README` prévus à cet effet, la totalité du TP doît être réalisé dans un notebook que vous joindrez à ce répository git en fin de séance.

### Données

#### Protéomiques

Les fiches UNIPROT des protéines étudiées dans la publication sont disponibles dans fichiers XML suivants:

- `data/Calderwood_EBV_proteome.xml`
- `data/Calderwood_Human_proteome.xml`

#### Interactomiques

Vous deverez recupérer les données d'interactions étudiées dans la publication grâce au protocole [PSICQUIC](https://psicquic.github.io/PsicquicSpec_1_4_Rest.html).

Ce protocole permet l'accès à distance à de nombreuses bases de données d'interactions protéine-protéine. Les interactions présentes dans ces bases de données sont obtenues par curation minutieuse de la littérature scientifique. Les interactions sont uniquement binaires (2 protéines) et toujours associées à la publication d'origine.
D'autres informations peuvent également être rapportées.
Le format standard de PSICQUIC, tabulé, est nommé [MITAB](https://psicquic.github.io/PSIMITAB.html).
Les requêtes sont formulées aux bases de données à l'aide d'un protocole REST obéissant à la syntaxe [MIQL](https://psicquic.github.io/MiqlReference.html)

Pour ce TP nous utiliserons la bases de données [Intact](https://www.ebi.ac.uk/intact/), dont le service PSICQUIC est accessible à l'URL: `http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search`.
Par exemple, les 100 premières interactions protéine-protéine humaines disponibles dans Intact sont accessibles via l'URL: `http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/species:human?firstResult=0&maxResults=100`

##### Quelles sont les significations des champs suivants du format MITAB 2.X?

Numero de champ | Signification Biologique|
 --- | --- 
1 | 
2 |
3 |
4 |
5 |
6 |

##### Utiliser le PMID de la publication pour récuperer les lignes MITAB des interactions rapportées dans l'étude.
Une librairie pratique pour manipuler des requêtes HTTP est [requests](https://requests.readthedocs.io/en/master/), eg:

```python
import requests
url = "https://mmsb.cnrs.fr/equipe/mobi/"

try:
    httpReq = requests.get(url, proxies=None)
except NameError:
    httpReq = requests.get(url)
ans = httpReq.text
```

##### Quelles techniques experimentales mesurent les interactions rapportées dans cette publication?

```

```

#### Extraction des deux sous-jeux d'interactions suivants:
- Interactions EBV-EBV
- Interactions EBV-Humaine


**Ne retenir que les lignes MITAB dans lesquelles chaque interactant possède un identifiant UNIPROT**

##### Proposer les expressions régulières et les champs auxquels les appliquer pour opérer les filtres suivants:
- Extraire les lignes MITAB impliquant uniquement des protéines d'EBV, quel est leur nombre ?
- Extraire les lignes MITAB impliquant des protéines humaines et des protéines d'EBV, quel est leur nombre ?

Jeux d'interactions | Numero de champ  | Expression(s) régulière(s) | Nombre d'interactions | Nombres d'interactants 
---                 |      ---         |          ---         |          ---          |         ---
**EBV-EBV**         |                  |      eg:`[\S]+`         |                       |
**EBV-Human**       |                  |      eg:`[\s]([\d]+)`   |                       |

##### Combien de protéines humaines et virales sont respectivement dans jeux d'interactions EBV-Human?

```

```

###### Pour la suite du travail assurez-vous d'avoir les deux jeux de données MITAB suivants

- MITAB EBV/EBV
- MITAB EBV/Human

Chacun ne contetant que des interactants référés par un numéro d'accession UNIPROT

### Construction du réseau d'interactions EBV/EBV

A l'aide des données MITAB et de la librarie [networkx](https://networkx.github.io/documentation/latest/
), représentez graphiquement un réseau dans lequel:

- les noeuds sont des identifiants UNIPROT

- les arêtes relient deux protéines en interaction

![Graphique](ebv_ebv_network_uniprot.png)

##### Décrivez brièvement ce réseau

```

```

Les noms de gènes sont parfois plus parlants que des accesseurs UNIPROT. A l'aide du fichier `./data/Calderwood_EBV_proteome.xml`  créez une table de conversion entre accesseur UNIPROT et nom de gène.

Pour vous aidez dans cette tâche, vous disposez du "parser" XML suivant qui étant donné un numéro d'accession UNIPROT et un document XML retourne un dictionaire d'informations concernant cette protéine. Le code permettant d'extraire l'information du nom de gène à été supprimée lors d'un `git push` malencontreux, à vous de le completer avant utilisation.

```python
from xml.etree.ElementTree import parse, dump, fromstring, register_namespace, ElementTree

# Utility functions
# Extracting All go terms relative to provided UNIPROT accessor
def goTerms(xmlEntry):
    ns = '{http://uniprot.org/uniprot}'
    goTerms = xmlEntry.findall(ns +'dbReference[@type="GO"]')
    goTermList = []
    for goT in goTerms:
        gID   = goT.attrib['id']
        gName = goT.find(ns +'property[@type="term"]').attrib['value']
        goTermList.append({"name" : gName, "ID" : gID})
    return goTermList

# Return information about provided UNIPROT accessor as python dictionary
def proteinDict(uniprotID, root):
    ns   = '{http://uniprot.org/uniprot}'

    data = { "accession" : uniprotID,
             "geneName" : None,
             "name" : None,
             "GOterms" : None
           }

    for entry in root.findall(ns+'entry'):
        accessions = entry.findall(ns+"accession")
        for acc in accessions:
            if acc.text == uniprotID: # entry is the node matching provided UNIPROT accessor
                e = entry.find(f"{ns}protein/{ns}recommendedName/{ns}fullName")
                if not e is None:
                    data["name"] = e.text
                e = "OUPSS##!!!"
                if not e is None:
                    data["geneName"] = e.text

                data["GOterms"] = goTerms(entry)
                return data
    raise ValueError(f"{uniprotID} nor found in XML document")
```

```python
# Test
tree = parse('./data/Calderwood_Human_proteome.xml')
root = tree.getroot()
proteinDict("Q53Y88", root)
```

Vous pouvez desormais dessiner le réseau dans lequel:

- les arêtes relient deux protéines en interaction
- les noeuds sont les noms des gènes correspondant aux protéines.

![Graphique](ebv_ebv_network_gene.png)

# Identification des cibles protéiques du virus.

Pour perturber et détourner à son profit le fonctionnement du lymphocyte, le matériel protéique du virus va interagir préferentiellement avec certains processus cellulaires.

En suivant, la méthode précédente dessiner le réseau dans lequel:

- les arêtes relient une protéine humaine et une protéine virale en interaction.
- les noeuds sont les noms des gènes correspondant aux protéines.
- les noeuds des protéines virales et humaines sont dessinés différement.

Vous pouvez jouer sur la taille de la figure et la constante de ressort *k* du rendu [spring_layout](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html) pour augmenter la lisibilité du graphique

![Graphique](ebv_human_network_gene.png)

###### Optionel

Les noeuds du réseau d'interaction peut aussi être divisés en deux partitions, humaine et virale. Vous disposez, ci-dessous d'un exemple de rendu graphique "bipartite" sur deux selections arbitraires de noeuds. Essayez de l'adapter au problème de représentation graphique du réseau d'interaction EBV/Humain

```python
import networkx as nx
plt.figure()
plt.axis('off')

G = nx.Graph()
G.add_edge("A", "C")
G.add_edge("B", "C")
X, Y = nx.bipartite.sets(G,  ["A", "B"])
pos = nx.bipartite_layout(G, X)
nx.draw_networkx_nodes(G, pos, node_color="blue", node_shape="o",nodelist=["A", "B"] )
nx.draw_networkx_nodes(G, pos, node_color="red", node_shape="o",nodelist=["C"] )
nx.draw_networkx_labels(G,pos,font_weight=800,font_color='black')
nx.draw_networkx_edges(G, pos, width=0.5)
```

##### Identification des processus biologiques ciblés par le virus

Nous allons observés les termes GO présents dans les protéines humaines les plus ciblées par EBV

1. Classer les noeuds, protéines humaines, par degré décroissant.
2. Lister les termes GO des protéines les plus connectées

- Quels termes reviennent fréquemment?

```

```

- Connaissant l'organisation de l'ontologie GO, quelles suggestions pourriez vous faire afin d'augmenter la qualité de l'analyse?

```

```

##### Bouquet final !

Bien souvent, un degré élevé d'interaction est requis pour visualiser des graphs complexes. Le JavaScript, notamment la librairie [D3](https://d3js.org/), est actuellement une solution technique de choix pour construire des représentations visuelles riches et interactives.

###### 1) Production d'un fichier JSON, encodant le graph d'interaction EBV/Human

- G2 est l'objet networkx *graph*
- `humanGeneLabels` et `ebvGeneLabels` sont des tables de conversions `{ clé[accesseur uniprot]:valeur[nom de gène] }`

```python
import json
nodes = [{'name': str(i),
          'id': humanGeneLabels[i] if i in humanGeneLabels else ebvGeneLabels[i],
          'type' : 'Human' if i in humanGeneLabels.keys() else 'EBV',
          'weight':G2.degree(i)
         }
         for i in G2.nodes()
        ]

links = []
for a,b in G2.edges() :
    _ = { "source" : None, "target" : None }

    for i, n in enumerate(nodes):
        if n['name'] == a:
            _["source"] = i
        elif n['name'] == b:
            _["target"] = i
    if  _["source"] is None or _["target"] is None:
        print(_)
        raise ValueError(a,b)
    links.append(_)
         
with open('graph.json', 'w') as f:
    json.dump({'nodes': nodes, 'links': links},
              f, indent=4)
```

###### 2) Construction d'une cellule de rendu Jupyter contenant la vue HTML
```python
%%html
<div id="d3-example"></div>
<style>
.node {stroke: #fff; stroke-width: 1.5px;}
.link {stroke: #999; stroke-opacity: .6;}
</style>
```

###### 3) Injection dans la cellule précedente du SVG représentant le réseau encodé dans `graph.json`

```javascript
%%javascript
// Fetch D3 library
const d3path = "https://d3js.org/d3.v6.min.js"
// Or load it locally, by default served from ~.jupyter/extensions
//const d3path = "d3.v6.min.js"
 console.log("Starting");
// We load the d3.js library from the Web.
require.config( { paths : { d3: d3path }
    } );

require(["d3"], function(d3) {
    console.log("Loading");
  // The code in this block is executed when the
  // d3.js library has been loaded.

  // First, we specify the size of the canvas
  // containing the visualization (size of the
  // <div> element).
  var width = 800, height = 800;

  // We create a color scale.
  var color = d3.scale.category10();

  // Create scale for node radius
  let rMin = 3, rMax = 12;
  let degMin = 1, degMax = 0;

  // We create a force-directed dynamic graph layout.
  var force = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);

  // In the <div> element, we create a <svg> graphic
  // that will contain our interactive visualization.
  var svg = d3.select("#d3-example").select("svg")
  if (svg.empty()) {
    svg = d3.select("#d3-example").append("svg")
          .attr("width", width)
          .attr("height", height);
  }

  // We load the JSON file.
  d3.json("graph.json", function(error, graph) {
    // In this block, the file has been loaded
    // and the 'graph' object contains our graph.
    console.log("Loading JSON")
    // We load the nodes and links in the
    // force-directed graph.
    console.dir(graph);
    force.nodes(graph.nodes);
    console.log("nodes loaded");
    force.links(graph.links)
    console.log("links loaded");
   
    // We create a <line> SVG element for each link
    // in the graph.
    var link = svg.selectAll(".link")
      .data(graph.links)
      .enter().append("line")
      .attr("class", "link");

    // We create a <circle> SVG element for each node
    // in the graph, and we specify a few attributes.
    var node = svg.selectAll(".node")
      .data(graph.nodes)
      .enter().append("circle")
      .attr("class", "node")
      //.attr("r", 5)  // radius
      .style("fill", function(d) {
          return  d.type === "EBV" ? "firebrick" : "steelblue"
      }).each(function(d) {
          degMax = d.weight > degMax ? d.weight : degMax;
      })
      .call(force.drag);
      // We parameterize the radius scale according to data
      let rScale = d3.scale.linear()
          .domain([degMin, degMax]) // unit: degree
          .range([rMin,rMax]); // unit: pixels
      node.attr("r", function(d) {
           return rScale(d.weight);
      });  // radius
      
    // The name of each node is the node number.
    node.append("title")
        .text(function(d) { return d.name; });

    // We bind the positions of the SVG elements
    // to the positions of the dynamic force-directed
    // graph, at each time step.
    force.on("tick", function() {
      link.attr("x1", function(d){return d.source.x})
          .attr("y1", function(d){return d.source.y})
          .attr("x2", function(d){return d.target.x})
          .attr("y2", function(d){return d.target.y});

      node.attr("cx", function(d){return d.x})
          .attr("cy", function(d){return d.y});
    });
      
    force.start();
    console.log("OK");
  });
},
function (err) {
    console.log("Could not load JS library " + err);
});
```
