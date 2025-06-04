import tools
import logging
import panel as pn
from panel.template import MaterialTemplate

logger = logging.getLogger(__name__)
logging.basicConfig(
    filename='dashboard.log',
    encoding='utf-8',
    level=logging.INFO
)

class AppState:
    def __init__(self):
        pass


def App() -> MaterialTemplate:
    # pn.extension('ipywidgets')
    app_state = AppState()

    template = pn.template.MaterialTemplate(
            title="title",
            header=[],
            main=[tools.visualize('10m', 'AL', 'SLOPE', '/media/volume/gabriel-geotiled/demo_test', downsample=10)], # put matplotlib figure here
            sidebar=[],
            header_background="#C4DEFC",
            busy_indicator=None,
            collapsed_sidebar=True
        )

    return template


try:
    app = App()
    app.servable()
except Exception as e:
    logger.error(f"dashboard could not be initialized {e}")