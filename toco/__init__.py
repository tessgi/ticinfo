from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import logging, os, __main__
logging.basicConfig()  # Avoid "No handlers could be found for logger" warning
logger = logging.getLogger(__name__)



from .toco import (toco)
PACKAGEDIR = os.path.dirname(os.path.abspath(__file__))