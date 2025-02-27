import numpy as np 
import os 
import pandas as pd

# Set range for the two variables of interest 
points = 5

# Diatoms 
ws1 = np.linspace(-1.38e-6, -1.38e-4, num=points)

# Microcystis 
ws2 = np.linspace(1.38e-5, 1.38e-3, num=points)

# Diatom growth rate
diatom = 1.388888888888889e-5
hab = 2.222222222222222e-6
pmax1 = np.logspace(0.25*diatom, 4*diatom, num=5)

# Microcystis growth rate
pmax2 = np.logspace(0.25*hab, 4*hab, num=5)

x1, x2, x3, x4  = np.meshgrid(ws1, ws2, pmax1, pmax2)

dataout = pd.DataFrame({'ws1': x1.flatten(), 'ws2': x2.flatten(), 'pmax1': x3.flatten(), 'pmax2': x4.flatten()})
dataout['fout_name'] = range(0, len(dataout))
output_csv = "parameter_space_ws_pmax.csv"
dataout.to_csv(output_csv, index=False)

# for p0 in pressure:

#     for ws0 in ws:
#         # print("Pressure = %f, ws = %f" % (p0, ws0))

#         data = {"pressure": [p0], "ws": [ws0]}
#         print(data)
#         df = pd.DataFrame(data)

#         # Append the DataFrame to the CSV file
#         file_exists = os.path.isfile(output_csv)
#         if file_exists:
#             mode = 'a'
#             header=False
#         else:
#             mode = 'w'
#             header=True

#         df.to_csv(output_csv, mode=mode, index=False, header=header)