import geopandas as gpd
from supply_demand import *
import configparser

config = configparser.ConfigParser(interpolation=None)
config.read('config.txt')

data_source_shortest = config['DATA_INPUT_TYPE']['TYPE']

data_output = config['DATA_OUTPUT_TYPE']['TYPE']

if data_source_shortest == 'FILE':
    shortest_gdf = gpd.read_file(config['OUTPUT_DATA_FILE']['KG_SHORTEST'])
    shortest_gdf['id_source'] = shortest_gdf['id_source'].astype(int)

elif data_source_shortest == 'DATABASE':
    conn = connect_database(config)

    schema_in = config['OUTPUT_DATA_DATABASE']['SCHEMA_OUT']
    table_shortest = config['OUTPUT_DATA_DATABASE']['TABLE_OUT_SHORTEST']

    kg_shortest_paths_query = f"SELECT * FROM {schema_in}.{table_shortest}"
    shortest_gdf = gpd.read_postgis(kg_shortest_paths_query, conn, geom_col='geometry')

# Assigns the res buildings to the kg.
kg_assigned = assign_shortest_paths(shortest_gdf)

if data_output == 'FILE':
    output = config['OUTPUT_DATA_FILE']['KG_RESULTS']
    joined_gdf = gpd.GeoDataFrame(shortest_gdf.merge(kg_assigned, on='s_t_ids', how='right')).to_file(output, driver='GeoJSON')

elif data_output == 'DATABASE':
    conn = connect_database(config)

    engine = create_engine('postgresql+psycopg2://', creator=lambda: conn)
    schema_out = config['OUTPUT_DATA_DATABASE']['SCHEMA_OUT']
    table_out = config['OUTPUT_DATA_DATABASE']['TABLE_OUT_RESULTS']

    joined_gdf = gpd.GeoDataFrame(shortest_gdf.merge(kg_assigned, on='s_t_ids', how='right'))

    joined_gdf.to_postgis(table_out, engine, schema=schema_out, if_exists='append', index=False, dtype={'geometry': 'geometry'})

