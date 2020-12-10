import bebi103

import iqplot

import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()

import holoviews as hv
hv.extension("bokeh")

bebi103.hv.set_defaults()



def plot_ecdf(df):
    p = iqplot.ecdf(
        df,
        q="time to catastrophe (s)",
        cats="concentration (μM)",
        style="staircase",
        conf_int=True,
    )

    p.legend.title = "concentration (μM)"
    p.legend.click_policy = "hide"
    bokeh.io.show(p)
    return p



def plot_stripbox(df):
    sb = iqplot.stripbox(
        data=df,
        q="time to catastrophe (s)",
        cats="concentration (μM)",
        jitter=True,
        top_level="box",
        marker_kwargs=dict(alpha=0.3)
    )
    bokeh.io.show(sb)
    return sb


