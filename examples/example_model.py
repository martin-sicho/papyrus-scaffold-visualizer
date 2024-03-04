"""
example_model.py

Created by: Martin Sicho
On: 22.05.23, 9:55
"""
from qsprpred import TargetTasks

from scaffviz.clustering.manifold import TSNE
from utils import fetch_example_models
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from scaffviz.depiction.plot import ModelPerformancePlot

# fetch two example models (single task classifier and single task regressor)
models, datasets = fetch_example_models(
    models=[
        RandomForestClassifier,
        RandomForestRegressor
    ],
    target_props=[
        [{
            "name": "pchembl_value_Median",
            "task": TargetTasks.SINGLECLASS,
            "th": [6.5]
        }], [{
            "name": "pchembl_value_Median",
            "task": TargetTasks.REGRESSION,
        }]
    ]
)
plot_types=(
    "errors", # plot mispredictions (difference between predictions and true labels)
    "splits", # plot folds and splits
    "predictions", # plot predicted values
    "labels", # plot original (true) labels
)
info = dict()
ports = [9000, 9100] # starting ports for each model
for plot_type in plot_types: # make a plot for each type (just an example, not recommended to do all at once)
    plot = ModelPerformancePlot(
        TSNE(), # use t-SNE for dimensionality reduction, does not recalculate if already done before on a data set
        models, # list of models to show the plot for
        datasets, # list of data sets used to fit the models
        ports,  # ports on localhost to use for each model performance plot (one for each model, must be unique)
        plot_type=plot_type, # type of the plot
        async_execution=True, # run the plot in a separate thread, set to False if you encounter problems
    )
    info_ = plot.make()
    info.update(info_)
    ports = (ports[0] + 1, ports[1] + 1)

# info about running plots
for port in info:
    print(f"The '{info[port]['plot_type']}' plot for model: '{info[port]['model'].name}' is running @ http://localhost:{port}")




