The content of the Readme.txt file is stuctured in the following way:
1. Description of neccessary python libraries
2. Description of the supply_demand algorythm files
3. Terms dictionary

1. The following packages need to be installed for the supply_demand algorythm:
- geopandas
- pandas
- networkx
- shapely
- psycopg2
- sqlalchemy
- configparser

2. There are four files for the algorythm to be utilized. Three python files and a configuration file are present: 
- supply_demand.py
- shortest_paths_demo.py
- assignment_demo.py
- config.txt - used to store credential data and parameter values

supply_demand.py
The first one contains the essential functions needed to implement the supply-demand analysis using a network and two types of point data representing accordingly the supply and the demand of the urban amenities.

shortest_paths_demo.py
The second one contains a procedure that implements the functions to create supply-demand (facility-building) pairs in a specified search radius.

assignment_demo.py
The third one contains the assignemnt procedure that assignes the reached demand buildings to the supply facilities taking into account the capacity of the facilities and their demand values as well as the fillment percentage.

config.txt
The fourth one is used to store credential data, to configure whether the data source and the data output is a file or a table in a database and to specify their names and location. This file is also used to set the parameter value for the Search radius.

3. The following terms can be found in the files
kg - kindergartens
res - residential