import pytest


def pytest_addoption(parser):
    """
    Add the --run-extra option to the pytest command line.
    """
    parser.addoption(
        '--run-extra',
        action='store_true',
        default=False,
        help='run extra tests',
    )


def pytest_configure(config):
    """
    Add the extra marker to the pytest configuration.
    """
    config.addinivalue_line('markers', 'extra: mark test as extra to run')


def pytest_collection_modifyitems(config, items):
    """
    Skip tests marked as extra if the --run-extra option is not provided.
    """
    if config.getoption('--run-extra'):
        return
    skip_extra = pytest.mark.skip(reason='need --run-extra option to run')
    for item in items:
        if 'extra' in item.keywords:
            item.add_marker(skip_extra)
