import logging

def setup_logging(basic_level="DEBUG", console_level="INFO", name=None):
    if name is None:
        name = __name__

    # Get the root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level=basic_level)
    # Setting matplotlib logger to warning
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("trimesh").setLevel(logging.ERROR)
    logging.getLogger("numba").setLevel(logging.ERROR)

    logger = logging.getLogger(name)
    logger.setLevel(level=basic_level)

    if root_logger.hasHandlers():
        root_logger.handlers.clear()

    # file handler for console logs
    all_file_handler = logging.FileHandler(
        "console.log",
        mode="a",
    )
    all_file_handler.setLevel(console_level)
    all_file_formatter = logging.Formatter(
        fmt="%(asctime)s-%(name)s-%(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    all_file_handler.setFormatter(all_file_formatter)

    # file handler specifically for debug logs
    debug_file_handler = logging.FileHandler(
        "report.log",
        mode="w",
    )
    debug_file_handler.setLevel(logging.DEBUG)
    debug_file_formatter = logging.Formatter(
        fmt="%(asctime)s-%(name)s-%(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    debug_file_handler.setFormatter(debug_file_formatter)

    # Setup a console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(console_level)

    # Add handlers to the root logger
    root_logger.addHandler(all_file_handler)
    root_logger.addHandler(debug_file_handler)
    root_logger.addHandler(console_handler)

    return logger

