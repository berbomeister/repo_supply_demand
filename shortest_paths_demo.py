import pandas as pd
import networkx as nx
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString
import psycopg2
from geoalchemy2 import Geometry, WKTElement
from sqlalchemy import *
import geopandas as gpd
from supply_demand import *
import configparser

# Start config file.
config = configparser.ConfigParser(interpolation=None)
config.read('config.txt')

data_source = config['DATA_INPUT_TYPE']['TYPE']

data_output_shortest = config['DATA_OUTPUT_TYPE']['TYPE']

if data_source == 'FILE':
    kg_gdf = gpd.read_file(config['INPUT_DATA_FILE']['KINDERGARTENS']).explode(index_parts=True) # kindergartens
    res_gdf = gpd.read_file(config['INPUT_DATA_FILE']['RESIDENTIAL_BUILDINGS']).explode(index_parts=True) # residential buildings
    ped_net_gdf = gpd.read_file(config['INPUT_DATA_FILE']['PED_NETWORK']).explode(index_parts=True) # pedestrian network


elif data_source == 'DATABASE':
    conn = connect_database(config)

    schema_in = config['INPUT_DATA_DATABASE']['SCHEMA_IN']
    table_kg = config['INPUT_DATA_DATABASE']['KINDERGARTENS_TABLE']
    table_res = config['INPUT_DATA_DATABASE']['RESIDENTIAL_BUILDINGS_TABLE']
    table_ped_net = config['INPUT_DATA_DATABASE']['PED_NETWORK_TABLE']

    kg_query = f"SELECT * FROM {schema_in}.{table_kg}"
    res_query = f"SELECT * FROM {schema_in}.{table_res}"
    ped_net_query = f"SELECT * FROM {schema_in}.{table_ped_net}"

    kg_gdf = gpd.read_postgis(kg_query, conn, geom_col='geometry').explode(index_parts=True)
    res_gdf = gpd.read_postgis(res_query, conn, geom_col='geometry').explode(index_parts=True)
    ped_net_gdf = gpd.read_postgis(ped_net_query, conn, geom_col='geometry').explode(index_parts=True)

search_radius = int(config['VARIABLES']['SEARCH_RADIUS'])  # meters

# Create graph from network dataset or database table.
G = create_graph(ped_net_gdf)

# Get kg nodes and res nodes from the buildings entrances datasets and slice them.
kg_nodes_buildings_sli_id = nodes_coordinates_buildings_sliced(kg_gdf, 'kg_id')
res_nodes_buildings_sli_id = nodes_coordinates_buildings_sliced(res_gdf, 'res_id')

# Get nodes from the graph.
nodes_graph_sli, nodes_graph_graphid = nodes_graph_sliced(G)

# Get kg nodes and res nodes id from the buildings entrances datasets and set it to kg and res nodes from the graph.
kg_graph_nodes_id = graph_nodes_match(kg_nodes_buildings_sli_id, nodes_graph_sli)
res_graph_nodes_id = graph_nodes_match(res_nodes_buildings_sli_id, nodes_graph_sli)

# Flip dictionary of id and coordinates for later use.
flipped_nodes_graph_graphid = flipped(nodes_graph_graphid)

# Match graphid of graph nodes with nodes id from buildings entrances datasets.
kg_graphid_id = graphid_id(kg_graph_nodes_id, nodes_graph_sli, nodes_graph_graphid)
res_graphid_id = graphid_id(res_graph_nodes_id, nodes_graph_sli, nodes_graph_graphid)

# Multi path shortest paths algo.
shortest_paths = shortest_path_algo(res_graph_nodes_id,kg_graphid_id,G,search_radius, nodes_graph_graphid)

# Export shortest paths data.
if data_output_shortest == 'FILE':
    output = config['OUTPUT_DATA_FILE']['KG_SHORTEST']
    sorted_shortest_paths, results_gdf = prepare_gdf_export(shortest_paths)
    export_shortest_paths(sorted_shortest_paths, results_gdf, G, kg_graphid_id, res_graphid_id, kg_gdf, res_gdf, flipped_nodes_graph_graphid, conn=None, config=None).to_file(output, driver='GeoJSON')

elif data_output_shortest == 'DATABASE':
    conn = connect_database(config)
    sorted_shortest_paths, results_gdf = prepare_gdf_export(shortest_paths)
    export_shortest_paths(sorted_shortest_paths, results_gdf, G, kg_graphid_id, res_graphid_id, kg_gdf, res_gdf, flipped_nodes_graph_graphid, conn, config)
