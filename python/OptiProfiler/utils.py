import logging


def get_logger(name=None, level=logging.INFO):
    logger = logging.getLogger(name)

    # Multiple calls to get_logger with the same name will return a reference
    # to the same logger. We do not create handlers if some already exist.
    if len(logger.handlers) == 0:
        logger.setLevel(level)

        # Attach a console handler (thread-safe).
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('[%(levelname)-8s] %(message)s'))
        logger.addHandler(handler)
    return logger
