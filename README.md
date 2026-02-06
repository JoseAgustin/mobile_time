# mobile_time

System to provide a temporal profiles for a gridded area by using  traffic activity data from 11,848 street segments, 8 vehicular types,  three different days and three road types.

## Information used
The following figure shows the data used for accomplish this task

![Input Information](/assets/images/diagrama.jpg "Input information")

The following tables describes the main information used in the code for compute the emissions and temporal profiles

   |  ID_time_period | Time period |
   |:---:|:---   |
   |  1 | Working Day |
   |  2 | Saturday |
   |  3 | Sunday   |
   | 11 | Week     |
   | 12 | Year     |

   |  type  | Description    | type | Description |
   |:---:    |:---            |:---: |:---          |
   | 11     | Automoviles     |  15  | Otros buses  |
   | 12     | Ligeros         |  16  | Medianos    |
   | 13     | Microbuses      |  17  | Pesados     |
   | 14     | Ruta 100        |  18  | Camiones Municipales |
   
 |  nef    | Description |
 |:---:    |:---             |
 | 1   | Automovil    |
 | 2   | Ligeros   |
 | 3  |  Microbuses   |
 | 4  |  Ruta 100, Pesados & Municipales TUV   |
 | 5   | Otros buses   |
 | 6  |  Medianos   |
     
  | kstype | Source classification | kstype |  Source classification |
  |:---:   |:---                     |:---: |:--- |
  |1 | Lateral              |    16  | Estación pesados             |
  |5 | Calle primaria       |    17  | Estación de autobuses        |
  |6 | Calle rápida         |    21  | Area residencial municipal   |
  |11 | area residencial    |    22  | Area res. e ind. municipal   |
  |12 | área res. e indus.  |    23  | Area pueblo municipal        |
  |13 | Area res. cerca centro|  24  | Estación de automóviles      |
  |14 | Area res. centro    |    25  | Area industrial Municipal    |
  |15 | Area pueblo         | | |

## Results

As a result the following products in a grid are generated:   daily emission  for NO, CO, SO2 and HC, hourly traffic temporal profile for the total fleet and for each vehicle type  and for each compound  for Mexico City for the 1990's.

The temporal profiles were obtains by computing  the ration between the hourly emission over the total emission per day per vehicular type and per chemical specie.

![CO emissions](/assets/images/COemis.gif "CO emissions")

![NO emissions](/assets/images/NOxemis.gif "NO emissions")
