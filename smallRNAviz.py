# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

# specify page formatting template
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# instantiate app object
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# specify app components, their data types and their starting values
app.layout = html.Div([
	dcc.Input(id="input-1", type="text", value="mouse"),
	dcc.Input(id="input-2", type="text", value="brain"),
	html.Div(id="output")
])

# specify which input values the app should listen for
@app.callback(
	Output("output", "children"),
	[Input("input-1", "value"), Input("input-2", "value")]
)

# the callback function (i.e. what happens after one of the inputs changes)
def update_output(input1, input2):
	return 'Input 1 is "{}" and Input 2 is "{}"'.format(input1, input2)

# launch the server when the program runs
if __name__ == "__main__":
	app.run_server(debug=True)

