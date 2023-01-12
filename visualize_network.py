import pygraphviz as pgv
import pandas as pd

edge_pallete = {'competition':"#fac631",
                'mutualism':"#0d0887"}

def build_network(edges, nodes, outfile, layout='dot'):
    edges = pd.read_csv(edges)
    nodes = pd.read_csv(nodes)
    g = pgv.AGraph(strict=False, directed=True)

    for node in nodes.itertuples():
        n = node.id
        l = node.label
        g.add_node(n,label=l)

    for edge in edges.itertuples():
        f= edge.source
        t= edge.target
        c= edge._4
        g.add_edge(f,t,color=edge_pallete[c])

    g.node_attr['shape']='box'
    g.node_attr['style']='filled'
    g.node_attr['arrowhead']='onormal'
    g.graph_attr['splines']='true'
    g.graph_attr['overlap']='false'
    g.graph_attr['size']='4,3!'
    g.graph_attr['dpi']='100'
    g.graph_attr['ratio']='fill'
    g.node_attr['fillcolor']='#ffffff'
    g.layout(layout)
    g.draw(outfile)

build_network('watches_edge_list.csv',
              'watches_node_list.csv',
              "figures/watches_graphviz.pdf",
              layout='circo')

build_network('law_edge_list.csv',
              'law_node_list.csv',
              "figures/law_graphviz.pdf")

build_network('realestate_edge_list.csv',
              'realestate_node_list.csv',
              "figures/realestate_graphviz.pdf",
              layout='circo')

build_network('data_edge_list.csv',
              'data_node_list.csv',
              "figures/data_graphviz.pdf")

build_network('redpill_edge_list.csv',
              'redpill_node_list.csv',
              "figures/redpill_graphviz.pdf")

build_network('mental_edge_list.csv',
              'mental_node_list.csv',
              "figures/mental_graphviz.pdf",
              layout='twopi')

build_network('cod_edge_list.csv',
              'cod_node_list.csv',
              "figures/cod_graphviz.pdf")
