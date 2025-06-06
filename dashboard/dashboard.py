import os
import yaml
import tools
import logging
import panel as pn
from matplotlib.figure import Figure
from panel.template import MaterialTemplate


# Initialize debug log file
logger = logging.getLogger(__name__)
logging.basicConfig(
    filename='dashboard.log',
    encoding='utf-8',
    level=logging.INFO
)


def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False  


def load_terrain_parameter(resolution, state, terrain_parameter, downsample, cmap, dpi, vmin, vmax):
    # Download the desired terrain parameter file if it does not yet exist
    tools.download(resolution, state, terrain_parameter)

    # Create figure for terrain parameter
    vmin = float(vmin) if is_float(vmin) else None
    vmax = float(vmax) if is_float(vmax) else None
    fig = tools.visualize(resolution, state, terrain_parameter, downsample, cmap, dpi, vmin, vmax)
    return fig


def App() -> MaterialTemplate:
    # Load in extension
    pn.extension()

    # Initialize pane variable
    pane = pn.pane.Matplotlib(Figure())
    
    # Load in configuration file
    config_path = os.getenv("CONFIG_PATH", f"./config.yaml")
    config = {}
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"could not initialize configuration, configuration path does not exists: {e}")
        return

    # Define widgets
    resolution = pn.widgets.Select(name='Resolution', options=list(config['resolutions'].keys()))
    state = pn.widgets.Select(name='State', options=list(config['states'].keys()))
    terrain_parameter = pn.widgets.Select(name='Terrain parameter', options=list(config['files'].keys()))
    downsample = pn.widgets.IntInput(value=1, start=1, end=100, name="Downsample Factor")
    dpi = pn.widgets.IntInput(value=300, start=100, end=500, name="DPI")
    cmap = pn.widgets.Select(name='Color Map', options=['viridis','plasma','inferno','magma','cividis'])
    vmin = pn.widgets.TextInput(value='', name='Normalized Minimum (Float)')
    vmax = pn.widgets.TextInput(value='', name='Normalized Maximum (Float)')
    button = pn.widgets.Button(name='Visualize', button_type='primary')

    # Bind different widgets to function to load and visualize terrain parameter
    load = pn.bind(load_terrain_parameter, resolution=resolution, state=state, terrain_parameter=terrain_parameter, downsample=downsample, cmap=cmap, dpi=dpi, vmin=vmin, vmax=vmax)

    # Create functionality for button
    def load_from_button(clicked):
        if clicked:
            fig = load()
            pane.object = fig
            
    load_button = pn.pane.Markdown(pn.bind(load_from_button, button, watch=True))
    
    # Organize controls in sidebar
    sidebar = pn.Column(pn.pane.Markdown('<h1>Options</h1> Larger data will take longer to download and visualize'), resolution, state, terrain_parameter, downsample, dpi, cmap, vmin, vmax, button, load_button)

    # Create template to display to dashboard
    template = pn.template.MaterialTemplate(
        title="GEOtiled Terrain Parameter Visualization Dashboard",
        header=[],
        main=[pane], # Panel (figure) goes here
        sidebar=[sidebar], # Controls goes here
        header_background="#203354",
        busy_indicator=None,
        collapsed_sidebar=True
    )

    return template


if __name__.startswith("bokeh"):
    try:
        app = App()
        app.servable()
    except Exception as e:
        logger.error(f"Dashboard could not be initialized {e}")