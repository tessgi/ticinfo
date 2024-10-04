from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import logging  # noqa: E401
import os

import __main__

PACKAGEDIR = os.path.dirname(os.path.abspath(__file__))
__version__ = "0.4.4"

import time
from threading import Event, Thread

from rich.console import Console
from rich.logging import RichHandler  # noqa: E402


def get_logger():
    """Configure and return a logger with RichHandler."""
    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.INFO)

    # # Add RichHandler
    # rich_handler = RichHandler(show_time=False, show_level=False, show_path=False, rich_tracebacks=True)
    # rich_handler.setFormatter(
    #     logging.Formatter("%(message)s")
    # )

    # logger.addHandler(rich_handler)
    #   return logger
    return TocoLogger("toco")


# Custom Logger with Rich
class TocoLogger(logging.Logger):
    def __init__(self, name, level=logging.INFO):
        super().__init__(name, level)
        console = Console()
        self.handler = RichHandler(
            show_time=False, show_level=False, show_path=False, console=console
        )
        self.handler.setFormatter(logging.Formatter("%(message)s"))
        self.addHandler(self.handler)
        self.spinner_thread = None
        self.spinner_event = None

    def start_spinner(self, message="Searching..."):
        if self.spinner_thread is None:
            self.spinner_event = Event()
            self.spinner_thread = Thread(target=self._spinner, args=(message,))
            self.spinner_thread.start()

    def stop_spinner(self):
        if self.spinner_thread is not None:
            self.spinner_event.set()
            self.spinner_thread.join()
            self.spinner_thread = None
            self.spinner_event = None

    def _spinner(self, message):
        with self.handler.console.status("[bold green]" + message) as status:
            while not self.spinner_event.is_set():
                time.sleep(0.1)


from .toco import toco  # noqa: E402, F401
