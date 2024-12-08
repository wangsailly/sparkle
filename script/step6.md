import csv
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import numpy as np
import random

###### #Initialize an empty list
data = []
node_link_path = "D:\MyCode\PythonCode\drawPic\edcs_genes.csv"

###### #Open the CSV file and read the contents
with open(node_link_path, newline='',encoding='utf-8') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:

        ###### ############
​        data.append(row)
data.pop(0) # Delete the column name

node_out = {} # The output occurrences of the node

###### #Count the occurrences
for row in data:
    if row[0] not in node_out:
        node_out[row[0]] = 1
    else:
        node_out[row[0]] += 1

###### #Count the nodes for drawing
node_list = []
for row in data:
    if node_out[row[0]] == 1:
        if "other-"+row[1] not in node_list:
           node_list.append("other-"+row[1])
    else:
        if row[0] not in node_list:
           node_list.append(row[0])
for row in data:
    if row[1] not in node_list:
        node_list.append(row[1])   

###### #Add HNSC node
if "HNSC" not in  node_list:
    node_list.append('HNSC') 
print(node_list)

node_list_index = {}
for idx, name in enumerate(node_list):
    node_list_index[name] = idx
source_list = []
taeget_list = []
value_list = []

###### #Add a connection
for row in data:
    if node_out[row[0]] == 1:
        source_list.append(node_list_index["other-"+row[1]])
        taeget_list.append(node_list_index[row[1]])
        value_list.append(1)
    else:
        source_list.append(node_list_index[row[0]])
        taeget_list.append(node_list_index[row[1]])
        value_list.append(1)
    ###### #Add HNSC connection
​    source_list.append(node_list_index[row[1]])
​    taeget_list.append(node_list_index['HNSC'])
​    value_list.append(1)

count_nodes = len(node_list)

###### #Use matplotlib to generate different colors based on individual preferences
cmap = plt.get_cmap('tab20')  
colors = [cmap(k / count_nodes) for k in range(count_nodes)]
colors = ['rgba' + str((int(r*255), int(g*255), int(b*255), a)) for r, g, b, a in colors]
random.shuffle(colors)

###### #####
links = {
    'source': source_list,  # Node index
    'target': taeget_list,  # Node index
    'value': value_list    
}

###### #Create a Sankey diagram
fig = go.Figure(data=[go.Sankey(
    node=dict(
      pad=50,
      thickness=50,
      line=dict(color="black", width=0.5),
      label=node_list,
      color=colors
    ),
    link=dict(
      source=links['source'],
      target=links['target'],
      value=links['value'],
  ))])

###### #Add node legend
for label, color in zip(node_list, colors):
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(
            size=10,
            color=color
        ),
        legendgroup=label,
        showlegend=True,
        name=label
    ))

fig.update_layout(title_text="Chemical-Gene-Disease Associations",width = 1080,height = 2000,font_size=10)
fig.show()
###### #save as PNG file
fig.write_image("Chemical-Gene-Disease Associations.png")

