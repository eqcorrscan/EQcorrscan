"""
Functions to handle logging for EQcorrscan

:copyright:
    EQcorrscan developers.

:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

import logging


class MockLoggingHandler(logging.Handler):
    """
    Mock logging handler to check for expected logs.

    Messages are available from an instance's ``messages`` dict, in order,
    indexed by a lowercase log level string (e.g., 'debug', 'info', etc.).

    Adapted from:
        https://stackoverflow.com/questions/899067/how-should-i-verify-a-
        log-message-when-testing-python-code-under-nose/20553331#20553331
    """

    def __init__(self, *args, **kwargs):
        self.messages = {'debug': [], 'info': [], 'warning': [], 'error': [],
                         'critical': []}
        super(MockLoggingHandler, self).__init__(*args, **kwargs)

    def emit(self, record):
        """
        Store a message from ``record`` in the instance's ``messages`` dict.
        """
        try:
            self.messages[record.levelname.lower()].append(record.getMessage())
        except Exception:
            self.handleError(record)

    def reset(self):
        self.acquire()
        try:
            for message_list in self.messages.values():
                message_list[:] = []
        finally:
            self.release()
