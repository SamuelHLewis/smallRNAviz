# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

import pandas as pd
import plotly.graph_objs as go

# specify page formatting template
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# read in data
sRNA_data = pd.DataFrame({
	"Length": range(17,30),
	"A_sense": [10,20,70,200,210,170,70,10,10,100,100,100,100],
	"C_sense": [10,20,70,200,210,170,70,10,10,100,100,100,100],
	"G_sense": [10,20,70,200,210,170,70,10,10,100,100,100,100],
	"U_sense": [10,20,70,200,210,170,70,10,1000,2000,2000,1000,100],
	"A_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-10,-100,-100,-100,-100],
	"C_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-10,-100,-100,-100,-100],
	"G_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-10,-100,-100,-100,-100],
	"U_antisense": [-10,-20,-70,-200,-210,-170,-70,-10,-1000,-2000,-2000,-1000,-100]
})

# instantiate app object
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# specify app components, their data types and their starting values
app.layout = html.Div([
	dcc.Graph(id = "graph-with-slider"),
	dcc.Slider(
		id = "length-slider",
		min = sRNA_data["Length"].min(),
		max = sRNA_data["Length"].max(),
		value = sRNA_data["Length"].min(),
		marks = {str(length): str(length) for length in sRNA_data["Length"].unique()}
	)
])

# specify which input values the app should listen for, and what it should output
@app.callback(
	Output("graph-with-slider", "figure"),
	[Input("length-slider", "value")]
)

# the callback function (i.e. what happens after one of the inputs changes)
def update_figure(selected_length):
	# create copy of data with desired length only
	filtered_data = sRNA_data[sRNA_data.Length == selected_length]
	# specify colour palette in one place, to make custom colours easier later on
	palette = ["green", "blue", "orange", "red"]
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

