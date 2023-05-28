"""
example_model.py

Created by: Martin Sicho
On: 22.05.23, 9:55
"""
from qsprpred.models.tasks import TargetTasks

from scaffviz.clustering.manifold import TSNE
from utils import fetch_example_models
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from scaffviz.depiction.plot import ModelPerformancePlot

models = fetch_example_models(
    models=[RandomForestClassifier, ExtraTreesClassifier],
    target_props=[{
        "name": "pchembl_value_Median",
        "task": TargetTasks.SINGLECLASS,
        "th": [6.5]
    }]
)
plot_types=(
    "errors",
    "splits",
    "predictions",
    "labels",
)
info = dict()
port_a = 9000
port_b = 9100
for plot_type in plot_types:
    plot = ModelPerformancePlot(TSNE(), models, plot_type=plot_type, async_execution=True, ports=(port_a, port_b))
    info_ = plot.make()
    info.update(info_)
    port_a += 1
    port_b += 1
print(info)




