# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ... import irsa
import numpy as np
import distutils.version as dv
import pytest

# this just wrong.  give up.
# @pytest.mark.skipif(dv.StrictVersion(np.__version__) <= dv.StrictVersion("1.4.1"))
# def test_trivial():
#     """ just make sure it doesn't raise anything
#     takes about 3-5 seconds"""
#     tbl = astroquery.irsa.query_gator_box('pt_src_cat','83.808 -5.391',300)
#
#     assert len(tbl) == 100 # at least, that's what I got...
#     return tbl
