# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
NEOCC Query Tool
================

Module to query the Near Earth Objects Coordination Centre (NEOCC).

"""

import os
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for 'ESANEOCC'
    """

    BASE_URL = 'https://' + os.getenv('NEOCC_PORTAL_IP',
                                      default='neo.ssa.esa.int')

    API_URL = _config.ConfigItem(BASE_URL +
                                 '/PSDB-portlet/download?file=')

    EPHEM_URL = _config.ConfigItem(BASE_URL +
                                   '/PSDB-portlet/ephemerides?des=')

    SUMMARY_URL = _config.ConfigItem(BASE_URL +
                                     '/search-for-asteroids?sum=1&des=')

    TIMEOUT = 60

    SSL_CERT_VERIFICATION = bool(int(os.getenv('SSL_CERT_VERIFICATION',
                                               default="1")))


conf = Conf()

from .core import neocc, ESAneoccClass

__all__ = ['neocc', 'ESAneoccClass', 'Conf', 'conf']
