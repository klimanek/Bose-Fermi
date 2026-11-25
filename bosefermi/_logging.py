"""Internal logging setup.

Each validator / matcher / analyzer instance owns the configuration of its
logger, so a ``debug=True`` call cannot leave a stream handler behind that
makes later ``debug=False`` calls noisy (see ROADMAP P1-4).
"""

import logging

_NULL_HANDLER = logging.NullHandler()


def make_logger(name: str, debug: bool, prefix: str) -> logging.Logger:
    """Return ``name`` configured for this instance.

    With ``debug`` the logger emits DEBUG records to stderr with a short
    ``prefix``; otherwise it is silenced through a ``NullHandler``.
    """
    logger = logging.getLogger(name)
    logger.handlers.clear()
    logger.propagate = False

    if debug:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(f"{prefix}: %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
    else:
        logger.addHandler(_NULL_HANDLER)
        logger.setLevel(logging.WARNING)

    return logger
