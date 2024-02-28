import pandas as pd

########################
#### Get the data frame
########################

# Convert the dictionary to a DataFrame
df = pd.DataFrame.from_dict(environments, orient='index')
# Environtment being the dictionary output from the pdb reading script.

# Save the DataFrame to a CSV file
df.to_csv('output.csv', index=True)

# Save the DataFrame to an Excel file
# df.to_excel('output.xlsx', index=True)
# Para esto hay que instalar openpyxl también, así que yo creo que con csv vamos bien.