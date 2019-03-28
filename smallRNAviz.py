# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import base64
import io
import pandas as pd
import plotly.graph_objs as go
import pysam
from Bio.Seq import Seq

# specify page formatting template
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# # read in data
# sRNA_data = pd.DataFrame({
# 	"Length": range(17,30),
# 	"A_sense": [10,20,70,200,210,170,70,10,10,100,100,100,100],
# 	"C_sense": [10,20,70,200,210,170,70,10,10,100,100,100,100],
# 	"G_sense": [10,20,70,200,210,170,70,10,10,100,100,100,100],
# 	"U_sense": [10,20,70,200,210,170,70,10,1000,2000,2000,1000,100],
# 	"A_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-10,-100,-100,-100,-100],
# 	"C_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-10,-100,-100,-100,-100],
# 	"G_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-10,-100,-100,-100,-100],
# 	"U_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-1000,-2000,-2000,-1000,-100]
# })

# create dataframe without any read counts to establish plotting area limits
sRNA_data = pd.DataFrame({
	"Length": range(17,30)
})

# instantiate app object
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# specify app components, their data types and their starting values
app.layout = html.Div([
	html.Div([
		dcc.Upload(
			id = "upload-data",
			children = html.Div(["Drag and drop or ", html.A("Select file")]),
			style = {
				"width": "100%",
				"height": "60px",
				"lineHeight": "60px",
				"borderWidth": "1px",
				"borderStyle": "dashed",
				"borderRadius": "5px",
				"textAlign": "center",
				"margin": "10px"
			},
			multiple = False
		)
	]),
	html.Div([
		html.Div([
			dcc.RadioItems(
			id = "strand-plotting",
			options = [{"label": i, "value": i} for i in ["Separate strands", "Combined strands"]],
			value = "Separate strands",
			labelStyle = {"display": "inline-block"}
			)
		], style={'width': '48%', 'display': 'inline-block'})
	]),
	dcc.Graph(id = "sRNA_graph"),
	dcc.RangeSlider(
		id = "length-slider",
		min = sRNA_data["Length"].min(),
		max = sRNA_data["Length"].max(),
		value = [sRNA_data["Length"].min(), sRNA_data["Length"].max()],
		marks = {str(length): str(length) for length in sRNA_data["Length"].unique()}
	),
])

# function to parse a user CSV file
def bam_to_sRNA_counts(contents):
	"""
	This function takes a bam file, and returns a pandas dataframe of counts for sRNAs of each length, with separate counts for each strand and 5' base
	"""
	# set up dataframe for sRNA data
	sRNA_data = pd.DataFrame(columns = ["Sequence", "Strand"])
	# read file contents
	#content_string = contents.split(",")[1]
	#decoded = base64.b64decode(content_string)
	#bam_contents = pysam.AlignmentFile(io.StringIO(decoded.decode('utf-8')), "rb")
	bam_contents = pysam.AlignmentFile(contents, "rb")
	sense_seqs = []
	antisense_seqs = []
	for line in bam_contents:
		# isolate the cigar string (6th field)
		line_contents = str(line).split("\t")
		cigar = line_contents[5]
		# keep only mapped reads
		if cigar.endswith("M"):
			# isolate their strand and sequence fields
			strand = line_contents[1]
			seq = line_contents[9]
			if line_contents[1] == "16":
				sRNA_data = sRNA_data.append({"Sequence": str(Seq(seq).reverse_complement()), "Strand": "-"}, ignore_index = True)
			else:
				sRNA_data = sRNA_data.append({"Sequence": seq, "Strand": "+"}, ignore_index = True)
	# add a length column
	sRNA_data["Length"] = sRNA_data.apply(lambda row: len(row["Sequence"]), axis = 1)
	# add a 5' base column
	sRNA_data["5base"] = sRNA_data.apply(lambda row: row["Sequence"][0], axis = 1)
	# create blank dataframe to hold counts of each base-strand combination for each length
	sRNA_counts = pd.DataFrame({
		"Length": list(range(17,36)),
		"A_sense": [0 for i in range(17,36)],
		"C_sense": [0 for i in range(17,36)],
		"G_sense": [0 for i in range(17,36)],
		"U_sense": [0 for i in range(17,36)],
		"A_antisense": [0 for i in range(17,36)],
		"C_antisense": [0 for i in range(17,36)],
		"G_antisense": [0 for i in range(17,36)],
		"U_antisense": [0 for i in range(17,36)]
	})
	# set the index as the sRNA lengths
	sRNA_counts = sRNA_counts.set_index("Length")
	# go through the bam_contents_df dataframe, updating the counts in the sRNA_counts dataframe
	# outer loop iterates over each length
	for length in sRNA_counts.index:
		# grab all sRNAs that have this length
		sRNAs_of_desired_length = sRNA_data[sRNA_data["Length"]==length]
		# check the strand of each sRNA
		for index, row in sRNAs_of_desired_length.iterrows():
			if row["Strand"] == "+" and row["5base"] == "A":
				sRNA_counts.at[length, "A_sense"] += 1
			elif row["Strand"] == "+" and row["5base"] == "C":
				sRNA_counts.at[length, "C_sense"] += 1
			elif row["Strand"] == "+" and row["5base"] == "G":
				sRNA_counts.at[length, "G_sense"] += 1
			elif row["Strand"] == "+" and row["5base"] == "T":
				sRNA_counts.at[length, "U_sense"] += 1
			elif row["Strand"] == "-" and row["5base"] == "A":
				sRNA_counts.at[length, "A_antisense"] -= 1
			elif row["Strand"] == "-" and row["5base"] == "C":
				sRNA_counts.at[length, "C_antisense"] -= 1
			elif row["Strand"] == "-" and row["5base"] == "G":
				sRNA_counts.at[length, "G_antisense"] -= 1
			elif row["Strand"] == "-" and row["5base"] == "T":
				sRNA_counts.at[length, "U_antisense"] -= 1
	# put Length back as a column
	sRNA_counts.reset_index(inplace = True)
	return sRNA_counts

# specify which input values the app should listen for, and what it should output
@app.callback(
	Output("sRNA_graph", "figure"),
	[Input("upload-data", "contents"),
	Input("length-slider", "value"),
	Input("strand-plotting", "value")]
)

# the callback function (i.e. what happens after one of the inputs changes)
def update_figure(user_data, length_range, strand_plotting):
	# specify colour palette in one place, to make custom colours easier later on
	palette = ["green", "blue", "orange", "red"]
	# read in user data
	sRNA_data = bam_to_sRNA_counts(user_data)
	# create copy of data with desired length only
	filtered_data = sRNA_data[(sRNA_data.Length >= length_range[0]) & (sRNA_data.Length <= length_range[1])]
	# if combined strand plotting is specified...
	if strand_plotting == "Combined strands":
		# make a new dataframe of summed sense and antisense read counts for each base
		filtered_data = pd.DataFrame({
			"Length": filtered_data.Length,
			"A": filtered_data.A_sense + abs(filtered_data.A_antisense),
			"C": filtered_data.C_sense + abs(filtered_data.C_antisense),
			"G": filtered_data.G_sense + abs(filtered_data.G_antisense),
			"U": filtered_data.U_sense + abs(filtered_data.U_antisense)
		})
		# specify 4 stacked bars per size value
		A = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.A,
			name='A',
			marker=dict(
				color=palette[0],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		C = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.C,
			name='C',
			marker=dict(
				color=palette[1],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		G = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.G,
			name='G',
			marker=dict(
				color=palette[2],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		U = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.U,
			name='U',
			marker=dict(
				color=palette[3],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		data = [U, G, C, A]
		layout = go.Layout(
			barmode='relative'
		)
		fig = go.Figure(data=data, layout=layout)
		return fig
	else:	
		# create each 5' base strand as a separate bar
		A_sense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.A_sense,
			name='A (sense)',
			marker=dict(
				color=palette[0],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		C_sense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.C_sense,
			name='C (sense)',
			marker=dict(
				color=palette[1],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		G_sense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.G_sense,
			name='G (sense)',
			marker=dict(
				color=palette[2],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		U_sense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.U_sense,
			name='U (sense)',
			marker=dict(
				color=palette[3],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		A_antisense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.A_antisense,
			name='A (antisense)',
			marker=dict(
				color=palette[0],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		C_antisense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.C_antisense,
			name='C (antisense)',
			marker=dict(
				color=palette[1],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		G_antisense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.G_antisense,
			name='G (antisense)',
			marker=dict(
				color=palette[2],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		U_antisense = go.Bar(
			x=filtered_data.Length,
			y=filtered_data.U_antisense,
			name='U (antisense)',
			marker=dict(
				color=palette[3],
				line=dict(
					color='black',
					width=1.5,
				)
			)
		)
		
		data = [U_sense, G_sense, C_sense, A_sense, U_antisense, G_antisense, C_antisense, A_antisense]
		layout = go.Layout(
			barmode='relative'
		)
		fig = go.Figure(data=data, layout=layout)
		return fig

# launch the server when the program runs
if __name__ == "__main__":
	app.run_server(debug=True)

