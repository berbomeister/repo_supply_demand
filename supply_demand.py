import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString
import psycopg2
from geoalchemy2 import Geometry, WKTElement
from sqlalchemy import *
import configparser

def connect_database(config):

    conn = psycopg2.connect(
        host=config['DATABASE_CREDENTIALS']['HOST'],
        database=config['DATABASE_CREDENTIALS']['DATABASE'],
        user=config['DATABASE_CREDENTIALS']['USER'],
        password=config['DATABASE_CREDENTIALS']['PASSWORD'].strip("''"),
        port=config['DATABASE_CREDENTIALS']['PORT']
    )
    return conn

def create_graph(network_gdf):
    # Create graph from network dataset or database table.

    network_gdf['geometry'] = [LineString(x.coords) for x in network_gdf.geometry]
    network_gdf['start_pos'] = [x.coords[0] for x in network_gdf.geometry]
    network_gdf['end_pos'] = [x.coords[-1] for x in network_gdf.geometry]

    s_points = network_gdf.start_pos.append(network_gdf.end_pos).reset_index(drop=True)
    s_points = s_points.drop_duplicates()

    df_points = pd.DataFrame(s_points, columns=['start_pos'])
    df_points['FNODE_'] = df_points.index
    network_gdf = pd.merge(network_gdf, df_points, on='start_pos', how='inner')

    df_points = pd.DataFrame(s_points, columns=['end_pos'])
    df_points['TNODE_'] = df_points.index
    network_gdf = pd.merge(network_gdf, df_points, on='end_pos', how='inner')

    df_points.columns = ['pos', 'osmid']
    df_points[['x', 'y']] = df_points['pos'].apply(pd.Series)
    df_node_xy = df_points.drop('pos', 1)

    graph = nx.Graph(crs=network_gdf.crs)

    for node, data in df_node_xy.T.to_dict().items():
        graph.add_node(node, **data)

    # Add edges to graph
    for i, row in network_gdf.iterrows():
        dict_row = row.to_dict()
        # if 'geometry' in dict_row: del dict_row['geometry']
        graph.add_edge(u_of_edge=dict_row['FNODE_'], v_of_edge=dict_row['TNODE_'], **dict_row)

    return graph


def nodes_coordinates_buildings_sliced(buildings_gdf, id_column):
    # Get kg nodes and res nodes from the buildings entrances datasets and slice them for matching with the nodes of the graph that are their representations.
    # This is necessary because of the slight change of the fractional part of the coordinates after the decimal point.

    nodes = {tuple(buildings_gdf.geometry.iloc[i].coords)[0]: buildings_gdf[id_column][i][0] for i in range(len(buildings_gdf))}
    nodes_dict_sli = {}
    for node in nodes:
        node_sli = (float(str(node[0])[:8]), float(str(node[1])[:9]))
        nodes_dict_sli[node_sli] = nodes[node]

    return nodes_dict_sli


def nodes_graph_sliced(graph):
    # Get nodes from the graph and slice them for matching with the buildings entrances datasets.

    nodes_dict_sli = {}
    nodes_graphid = {}

    for k, v in graph.nodes(data=True):
        node_sli = (float(str(v['x'])[:8]), float(str(v['y'])[:9]))
        nodes_dict_sli[node_sli] = (v['x'], v['y'])
        nodes_graphid[(v['x'], v['y'])] = k

    return nodes_dict_sli, nodes_graphid


def graph_nodes_match(nodes_buildings_sli, nodes_graph_sli):
    # Get kg nodes and res nodes id from the buildings entrances datasets and set it to kg and res nodes from the graph.
    # This way the graph nodes receive buildings attributes.

    keys = set(nodes_buildings_sli.keys()).intersection(set(nodes_graph_sli.keys()))
    graph_nodes_match = {nodes_graph_sli[k]: nodes_buildings_sli[k] for k in keys}  # coord : id

    return graph_nodes_match


def flipped(dict):
    # Flip dictionary of id and coordinates.
    result = {value: key for key, value in dict.items()}
    return result


def graphid_id(graph_nodes, nodes_graph_sli, nodes_graph_graphid):
    # Match graphid of graph nodes with nodes id from buildings entrances datasets.

    keys = set(graph_nodes.keys()).intersection(set(nodes_graph_sli.values()))
    result = {nodes_graph_graphid[k]: graph_nodes[k] for k in keys}

    return result


def shortest_path_algo(res_graph_nodes_id, kg_graphid_id, graph, search_radius, nodes_graph_graphid):
    # Multi path shortest paths algo.

    shortest_paths = {}
    node_count = 0
    errors = []
    target_points = MultiPoint(list(res_graph_nodes_id.keys()))

    for s_node in kg_graphid_id:
        try:
            path = nx.single_source_dijkstra_path(graph, s_node, cutoff=search_radius, weight='length_m') # Check the correct field for the weight
        except nx.NetworkXNoPath:
            errors.append(list(s_node))

        s_node_point = Point(graph.nodes[s_node]['x'], graph.nodes[s_node]['y'])
        near_target_points = target_points.intersection(s_node_point.buffer(search_radius))
        near_target_nodes = [(point.x, point.y) for point in near_target_points]
        near_target_nodes_graphid = [v for k, v in nodes_graph_graphid.items() if k in near_target_nodes]

        target_nodes_buffer_final = list(set(path.keys()).intersection(set(near_target_nodes_graphid)))

        for t_node in target_nodes_buffer_final:
            a = s_node, t_node
            shortest_paths[a] = path[t_node]

        node_count += 1
        print(f"{node_count} of {len(kg_graphid_id)} |", end=" ")

    return shortest_paths

def prepare_gdf_export(shortest_paths):
    sorted_shortest_paths = dict(sorted(shortest_paths.items()))
    results_gdf = gpd.GeoDataFrame(
        columns=['geometry', 'id', 'length_m', 'id_source', 'id_target', 'kids_kg_ca',
                 'kids_res_t', 'graphid_so', 'graphid_ta', 'coor_so',
                 'coor_ta'],
        dtype='object')

    return sorted_shortest_paths, results_gdf

def export_shortest_paths(sorted_shortest_paths, results_gdf, graph, kg_graphid_id, res_graphid_id, kg_gdf, res_gdf, flipped_nodes_graph_graphid, conn=None, config=None):
    # Loop through sorted_shortest_paths.
    count = 0
    errors = []
    prev_kg = next(iter(sorted_shortest_paths), None)[0]

    if conn is not None:
        engine = create_engine('postgresql+psycopg2://', creator=lambda: conn)
        schema_out = config['OUTPUT_DATA_DATABASE']['SCHEMA_OUT']
        table_out = config['OUTPUT_DATA_DATABASE']['TABLE_OUT_SHORTEST']

    for tup in sorted_shortest_paths:
        current_kg_node_graph = tup[0]

        if current_kg_node_graph != prev_kg and conn is not None:
            results_gdf["geometry"] = [MultiLineString([feature]) if isinstance(feature, LineString) else feature for feature in results_gdf["geometry"]]
            results_gdf.to_postgis(table_out, engine, schema=schema_out, if_exists='append', index=False)
            results_gdf = results_gdf[0:0]

        edges_in_shortest_path = []

        for i in range(len(sorted_shortest_paths[tup])-1):
            edges_in_shortest_path.append(graph.get_edge_data(sorted_shortest_paths[tup][i], sorted_shortest_paths[tup][i+1]))

        try:
            shortest_gdf = gpd.GeoDataFrame(edges_in_shortest_path, geometry='geometry', crs=7801)
        except ValueError:
            errors.append(list(current_kg_node_graph))

        shortest_gdf = gpd.GeoDataFrame.dissolve(shortest_gdf, aggfunc='sum')
        shortest_gdf['id_source'] = kg_graphid_id[tup[0]]
        shortest_gdf['id_target'] = res_graphid_id[tup[1]]
        shortest_gdf['s_t_ids'] = str(kg_graphid_id[tup[0]]) + "_" + str(int(res_graphid_id[tup[1]]))

        shortest_gdf['kids_kg_ca'] = float(kg_gdf.loc[kg_gdf['kg_id'] == kg_graphid_id[tup[0]], 'kids_all'])
        shortest_gdf['kids_res_t'] = float(res_gdf.loc[res_gdf['res_id'] == res_graphid_id[tup[1]], 'kids_count'])

        shortest_gdf['graphid_so'] = str(tup[0])
        shortest_gdf['graphid_ta'] = str(tup[1])
        shortest_gdf['coor_so'] = str(flipped_nodes_graph_graphid[tup[0]])
        shortest_gdf['coor_ta'] = str(flipped_nodes_graph_graphid[tup[1]])
        results_gdf = pd.concat([results_gdf, shortest_gdf], ignore_index=True)

        prev_kg = tup[0]

        count += 1
        print(f"{count} of {len(sorted_shortest_paths)} id = {kg_graphid_id[tup[0]]} |", end=" ")

    if conn is not None:
        results_gdf["geometry"] = [MultiLineString([feature]) if isinstance(feature, LineString) else feature for feature in results_gdf["geometry"]]
        results_gdf.to_postgis(table_out, engine, schema=schema_out, if_exists='append', index=False, dtype={'geometry': 'geometry'})

    return results_gdf

def assign_shortest_paths(shortest_gdf):
    # Assigns the res buildings to the kg.
    sorted_shortest_gdf = shortest_gdf.sort_values(['id_source', 'length_m'], ascending=[True, True])

    assignment_df = pd.DataFrame(columns=['id_target', 'assigned', 's_t_ids'])
    source_capacity = 0
    source_count = 0

    for source in sorted_shortest_gdf['id_source'].unique():
        source_capacity = sorted_shortest_gdf[sorted_shortest_gdf['id_source'] == source].reset_index().loc[0]['kids_kg_ca']
        for index, row in sorted_shortest_gdf[sorted_shortest_gdf['id_source'] == source].iterrows():
            target = int(row['id_target'])
            if target in assignment_df['id_target'].values:
                continue
            elif row['kids_res_t'] > source_capacity:
                continue
            else:
                assignment_df = assignment_df.append({
                    # 'id_source': int(source),
                    'id_target': int(target),
                    'assigned': True,
                    # 'length_m': row['length_m'],
                    # 'kids_res_t': row['kids_res_t'],
                    's_t_ids': str(source) + "_" + str(int(target)),
                    # 'geometry': row['geometry']
                }, ignore_index=True)
                source_capacity -= row['kids_res_t']

        source_count += 1
        print(f"{source_count} out of {len(sorted_shortest_gdf['id_source'].unique())} | ")

    assignment_df = assignment_df.drop(columns=['id_target'])

    return assignment_df


