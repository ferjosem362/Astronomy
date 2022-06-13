# Licensed under a 3-clause BSD style license - see LICENSE.rst

from cmath import exp
from collections import OrderedDict
import os
import requests
import pytest
import tempfile
import textwrap
import urllib

import astropy.coordinates as coord
from astropy.io import fits
import astropy.io.votable as votable
import astropy.units as u
from astropy.table import Table
import astropy.utils.data as aud
from astropy.utils.exceptions import AstropyDeprecationWarning

import unittest
from ...utils import chunk_read, chunk_report, class_or_instance, commons
from ...utils.process_asyncs import async_to_sync_docstr, async_to_sync
from ...utils.docstr_chompers import remove_sections, prepend_docstr_nosections


class SimpleQueryClass:

    @class_or_instance
    def query(self):
        """ docstring """
        if self is SimpleQueryClass:
            print("Calling query as class method")
            return "class"
        else:
            print("Calling query as instance method")
            return "instance"


@pytest.mark.remote_data
def test_chunk_read():
    datasize = 50000
    response = urllib.request.urlopen(f'http://httpbin.org/stream-bytes/{datasize}')
    C = chunk_read(response, report_hook=chunk_report)
    assert len(C) == datasize


def test_class_or_instance():
    assert SimpleQueryClass.query() == "class"
    U = SimpleQueryClass()
    assert U.query() == "instance"
    assert SimpleQueryClass.query.__doc__ == " docstring "


@pytest.mark.parametrize(('coordinates'),
                         [coord.SkyCoord(ra=148, dec=69, unit=(u.deg, u.deg)),
                          ]
                         )
def test_parse_coordinates_1(coordinates):
    c = commons.parse_coordinates(coordinates)
    assert c is not None


@pytest.mark.remote_data
@pytest.mark.parametrize(('coordinates'),
                         ["00h42m44.3s +41d16m9s",
                          "m81"])
def test_parse_coordinates_2(coordinates):
    c = commons.parse_coordinates(coordinates)
    assert c is not None


def test_parse_coordinates_3():
    with pytest.raises(Exception):
        commons.parse_coordinates(9.8 * u.kg)


def test_parse_coordinates_4():
    # Regression test for #1251
    coordinates = "251.51 32.36"
    c = commons.parse_coordinates(coordinates)
    assert c.to_string() == coordinates


def test_send_request_post(monkeypatch):
    def mock_post(url, data, timeout, headers={}, status_code=200):
        class SpecialMockResponse:

            def __init__(self, url, data, headers, status_code):
                self.url = url
                self.data = data
                self.headers = headers
                self.status_code = status_code

            def raise_for_status(self):
                pass

        return SpecialMockResponse(url, data, headers=headers,
                                   status_code=status_code)
    monkeypatch.setattr(requests, 'post', mock_post)

    with pytest.warns(AstropyDeprecationWarning):
        response = commons.send_request('https://github.com/astropy/astroquery',
                                            data=dict(msg='ok'), timeout=30)
    assert response.url == 'https://github.com/astropy/astroquery'
    assert response.data == dict(msg='ok')
    assert 'astroquery' in response.headers['User-Agent']
    assert response.headers['User-Agent'].endswith("_testrun")


def test_send_request_get(monkeypatch):
    def mock_get(url, params, timeout, headers={}, status_code=200):
        req = requests.Request(
            'GET', url, params=params, headers=headers).prepare()
        req.status_code = status_code
        req.raise_for_status = lambda: None
        return req
    monkeypatch.setattr(requests, 'get', mock_get)
    with pytest.warns(AstropyDeprecationWarning):
        response = commons.send_request('https://github.com/astropy/astroquery',
                                        dict(a='b'), 60, request_type='GET')
    assert response.url == 'https://github.com/astropy/astroquery?a=b'


def test_quantity_timeout(monkeypatch):
    def mock_get(url, params, timeout, headers={}, status_code=200):
        req = requests.Request(
            'GET', url, params=params, headers=headers).prepare()
        req.status_code = status_code
        req.raise_for_status = lambda: None
        return req
    monkeypatch.setattr(requests, 'get', mock_get)
    with pytest.warns(AstropyDeprecationWarning):
        response = commons.send_request('https://github.com/astropy/astroquery',
                                        dict(a='b'), 1 * u.min,
                                        request_type='GET')
    assert response.url == 'https://github.com/astropy/astroquery?a=b'


def test_send_request_err():
    with pytest.raises(ValueError):
        with pytest.warns(AstropyDeprecationWarning):
            commons.send_request('https://github.com/astropy/astroquery',
                                 dict(a='b'), 60, request_type='PUT')


col_1 = [1, 2, 3]
col_2 = [0, 1, 0, 1, 0, 1]
col_3 = ['v', 'w', 'x', 'y', 'z']
# table t1 with 1 row and 3 cols
t1 = Table([col_1[:1], col_2[:1], col_3[:1]], meta={'name': 't1'})
# table t2 with 3 rows and 1 col
t2 = Table([col_1], meta={'name': 't2'})
# table t3 with 3 cols and 3 rows
t3 = Table([col_1, col_2[:3], col_3[:3]], meta={'name': 't3'})


def test_TableDict():
    in_list = create_in_odict([t1, t2, t3])
    table_list = commons.TableList(in_list)
    repr_str = table_list.__repr__()
    assert repr_str == ("TableList with 3 tables:"
                        "\n\t'0:t1' with 3 column(s) and 1 row(s) "
                        "\n\t'1:t2' with 1 column(s) and 3 row(s) "
                        "\n\t'2:t3' with 3 column(s) and 3 row(s) ")


def test_TableDict_print_table_list(capsys):
    in_list = create_in_odict([t1, t2, t3])
    table_list = commons.TableList(in_list)
    table_list.print_table_list()
    out, err = capsys.readouterr()
    assert out == ("TableList with 3 tables:"
                   "\n\t'0:t1' with 3 column(s) and 1 row(s) "
                   "\n\t'1:t2' with 1 column(s) and 3 row(s) "
                   "\n\t'2:t3' with 3 column(s) and 3 row(s) \n")


def create_in_odict(t_list):
    return OrderedDict([(t.meta['name'], t) for t in t_list])


def test_suppress_vo_warnings(recwarn):
    commons.suppress_vo_warnings()
    votable.exceptions.warn_or_raise(votable.exceptions.W01)
    votable.exceptions.warn_or_raise(votable.exceptions.VOTableChangeWarning)
    votable.exceptions.warn_or_raise(votable.exceptions.VOWarning)
    votable.exceptions.warn_or_raise(votable.exceptions.VOTableSpecWarning)
    with pytest.raises(Exception):
        recwarn.pop(votable.exceptions.VOWarning)


docstr1 = """
        Query the Vizier service for a specific catalog

        Parameters
        ----------
        catalog : str or list, optional
            The catalog(s) that will be retrieved

        Returns
        -------
        response : `~request.response`
            Returned if asynchronous method used
        result : `~astroquery.utils.common.TableList`
            The results in a list of `astropy.table.Table`.
        """

docstr1_out = textwrap.dedent("""
        Queries the service and returns a table object.

        Query the Vizier service for a specific catalog

        Parameters
        ----------
        catalog : str or list, optional
            The catalog(s) that will be retrieved

        Returns
        -------
        table : A `~astropy.table.Table` object.
        """)

docstr2 = """
        Search Vizier for catalogs based on a set of keywords, e.g. author name

        Parameters
        ----------
        keywords : list or string
            List of keywords, or space-separated set of keywords.
            From `Vizier <http://vizier.u-strasbg.fr/doc/asu-summary.htx>`_:
            "names or words of title of catalog. The words are and'ed, i.e.
            only the catalogues characterized by all the words are selected."

        Returns
        -------
        Dictionary of the "Resource" name and the VOTable resource object.
        "Resources" are generally publications; one publication may contain
        many tables.

        Examples
        --------
        >>> from astroquery.vizier import Vizier
        >>> catalog_list = Vizier.find_catalogs('Kang W51')
        >>> print(catalog_list)
        {u'J/ApJ/706/83': <astropy.io.votable.tree.Resource at 0x108d4d490>,
         u'J/ApJS/191/232': <astropy.io.votable.tree.Resource at 0x108d50490>}
        >>> print({k:v.description for k,v in catalog_list.items()})
        {u'J/ApJ/706/83': u'Embedded YSO candidates in W51 (Kang+, 2009)',
         u'J/ApJS/191/232': u'CO survey of W51 molecular cloud (Bieging+, 2010)'}
        """

# note that the blank line between "Queries..." and "Search..." is necessary
# for sphinx parsing of the docs
docstr2_out = textwrap.dedent("""
        Queries the service and returns a dict object.

        Search Vizier for catalogs based on a set of keywords, e.g. author name

        Parameters
        ----------
        keywords : list or string
            List of keywords, or space-separated set of keywords.
            From `Vizier <http://vizier.u-strasbg.fr/doc/asu-summary.htx>`_:
            "names or words of title of catalog. The words are and'ed, i.e.
            only the catalogues characterized by all the words are selected."

        Examples
        --------
        >>> from astroquery.vizier import Vizier
        >>> catalog_list = Vizier.find_catalogs('Kang W51')
        >>> print(catalog_list)
        {u'J/ApJ/706/83': <astropy.io.votable.tree.Resource at 0x108d4d490>,
         u'J/ApJS/191/232': <astropy.io.votable.tree.Resource at 0x108d50490>}
        >>> print({k:v.description for k,v in catalog_list.items()})
        {u'J/ApJ/706/83': u'Embedded YSO candidates in W51 (Kang+, 2009)',
         u'J/ApJS/191/232': u'CO survey of W51 molecular cloud (Bieging+, 2010)'}

        Returns
        -------
        dict : A `dict` object.
        """)


def test_process_async_docs():
    assert async_to_sync_docstr(docstr1) == docstr1_out
    assert async_to_sync_docstr(docstr2, returntype='dict') == docstr2_out


class Dummy:

    def do_nothing_async(self):
        """ docstr """
        pass


def test_async_to_sync(cls=Dummy):
    newcls = async_to_sync(Dummy)
    assert hasattr(newcls, "do_nothing")


docstr3 = """
    Parameters
    ----------
    first_param

    Returns
    -------
    Nothing!

    Examples
    --------
    Nada
"""

docstr3_out = """
    Examples
    --------
    Nada
"""


def test_return_chomper(doc=docstr3, out=docstr3_out):
    assert (remove_sections(doc, sections=['Returns', 'Parameters']) ==
            [x.lstrip() for x in out.split('\n')])


def dummyfunc1():
    """
    Returns
    -------
    Nothing!

    Examples
    --------
    Nada
    """
    pass


def dummyfunc2():
    """
    Returns
    -------
    Nothing!
    """
    pass


docstr4 = """
    Blah Blah Blah

    Returns
    -------
    nothing

    Examples
    --------
    no_examples_at_all
"""

docstr4_out1 = """
    Blah Blah Blah

    Returns
    -------
    Nothing!

    Examples
    --------
    Nada
"""

docstr4_out2 = """
    Blah Blah Blah

    Returns
    -------
    Nothing!
"""


@pytest.mark.parametrize("func, out", [(dummyfunc1, docstr4_out1),
                                       (dummyfunc2, docstr4_out2)])
def test_prepend_docstr(func, out, doc=docstr4):
    fn = prepend_docstr_nosections(doc, sections=['Returns', 'Examples'])(func)
    assert fn.__doc__ == textwrap.dedent(out)


@async_to_sync
class DummyQuery:

    @class_or_instance
    def query_async(self, *args, **kwargs):
        """ docstr"""
        if kwargs['get_query_payload']:
            return dict(msg='payload returned')
        return 'needs to be parsed'

    @class_or_instance
    def _parse_result(self, response, verbose=False):
        return response


def test_payload_return(cls=DummyQuery):
    result = DummyQuery.query(get_query_payload=True)
    assert isinstance(result, dict)
    result = DummyQuery.query(get_query_payload=False)
    assert isinstance(result, str)


fitsfilepath = os.path.join(os.path.dirname(__file__),
                            '../../sdss/tests/data/emptyfile.fits')


@pytest.fixture
def patch_getreadablefileobj(request):
    # Monkeypatch hack: ALWAYS treat as a URL
    _is_url = aud._is_url
    aud._is_url = lambda x: True

    if not commons.ASTROPY_LT_4_3:
        _try_url_open = aud._try_url_open
        aud._try_url_open = lambda x, **kwargs: MockRemote(x, **kwargs)

    _urlopen = urllib.request.urlopen
    _urlopener = urllib.request.build_opener
    _urlrequest = urllib.request.Request
    filesize = os.path.getsize(fitsfilepath)

    class MockRemote:
        def __init__(self, fn, *args, **kwargs):
            self.file = open(fn, 'rb')

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            self.close()

        def info(self):
            return {'Content-Length': filesize}

        def read(self, *args):
            return self.file.read(*args)

        def close(self):
            self.file.close()

    def monkey_urlopen(x, *args, **kwargs):
        print("Monkeyed URLopen")
        return MockRemote(fitsfilepath, *args, **kwargs)

    def monkey_builder(tlscontext=None):
        mock_opener = type('MockOpener', (object,), {})()
        mock_opener.open = lambda x, **kwargs: MockRemote(fitsfilepath, **kwargs)
        return mock_opener

    def monkey_urlrequest(x, *args, **kwargs):
        # urlrequest allows passing headers; this will just return the URL
        # because we're ignoring headers during mocked actions
        print("Monkeyed URLrequest")
        return x

    aud.urllib.request.Request = monkey_urlrequest
    aud.urllib.request.urlopen = monkey_urlopen
    aud.urllib.request.build_opener = monkey_builder
    urllib.request.urlopen = monkey_urlopen
    urllib.request.build_opener = monkey_builder

    def closing():
        aud._is_url = _is_url

        if not commons.ASTROPY_LT_4_3:
            aud._try_url_open = _try_url_open

        urllib.request.urlopen = _urlopen
        aud.urllib.request.urlopen = _urlopen
        urllib.request.build_opener = _urlopener
        aud.urllib.request.build_opener = _urlopener
        aud.urllib.request.Request = _urlrequest

    request.addfinalizer(closing)


def test_filecontainer_save(patch_getreadablefileobj):
    ffile = commons.FileContainer(fitsfilepath, encoding='binary')
    temp_dir = tempfile.mkdtemp()
    empty_temp_file = f"{temp_dir}{os.sep}test_emptyfile.fits"
    ffile.save_fits(empty_temp_file)
    assert os.path.exists(empty_temp_file)


def test_filecontainer_get(patch_getreadablefileobj):
    ffile = commons.FileContainer(fitsfilepath, encoding='binary')
    ff = ffile.get_fits()
    assert isinstance(ff, fits.HDUList)


@pytest.mark.parametrize(('coordinates', 'expected'),
                         [("5h0m0s 0d0m0s", True),
                          ("m1", False)
                          ])
def test_is_coordinate(coordinates, expected):
    out = commons._is_coordinate(coordinates)
    assert out == expected


@pytest.mark.parametrize(('radius'),
                         [0.01*u.deg, '0.01 deg', 0.01*u.arcmin]
                         )
def test_radius_to_unit(radius):
    c = commons.radius_to_unit(radius)
    assert c is not None


# This method will be used by the mock in test_get_access_url to replace requests.get
def _mocked_requests_get(*args, **kwargs):
    class MockResponse:
        def __init__(self, text, status_code):
            self.text = text
            self.status_code = status_code

        def text(self):
            return self.text

        def raise_for_status(self):
            return None

    example_1_capabilities_document = """<?xml version="1.0" encoding="UTF-8"?>
<vosi:capabilities xmlns:vosi="http://www.ivoa.net/xml/VOSICapabilities/v1.0" xmlns:vs="http://www.ivoa.net/xml/VODataService/v1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<capability standardID="ivo://ivoa.net/auth#example">
<interface xsi:type="vs:ParamHTTP">
    <accessURL use="full">https://example.com/service1/endpoint</accessURL>
</interface>
</capability>
</vosi:capabilities>"""

    example_2_capabilities_document = """<?xml version="1.0" encoding="UTF-8"?>
<vosi:capabilities xmlns:vosi="http://www.ivoa.net/xml/VOSICapabilities/v1.0" xmlns:vs="http://www.ivoa.net/xml/VODataService/v1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<capability standardID="ivo://ivoa.net/auth#example">
<interface xsi:type="vs:ParamHTTP">
    <accessURL use="full">https://example.com/service2/endpoint</accessURL>
</interface>
</capability>
</vosi:capabilities>"""

    if args[0] == 'https://example.com/reg/resource-caps':
        return MockResponse("ivo://example.com/service1 = https://example.com/service1/capabilities\nivo://example.com/service2 = https://example.com/service2/capabilities", 200)
    elif args[0] == 'https://example.com/service1/capabilities':
        return MockResponse(example_1_capabilities_document, 200)
    elif args[0] == 'https://example.com/service2/capabilities':
        return MockResponse(example_2_capabilities_document, 200)

    pytest.fail('Should not get here.')


@pytest.mark.parametrize('service, capability, expected',
                         [
                            ('ivo://example.com/service1', 'ivo://ivoa.net/auth#example', 'https://example.com/service1/endpoint'),
                            ('ivo://example.com/service2', 'ivo://ivoa.net/auth#example', 'https://example.com/service2/endpoint'),
                            ('bogus', '', ''),
                            ('', '', '')
                         ])
def test_get_access_url(service, capability, expected):
    """
    Test for get_access_url to obtain a download URL from a Capabilities document
    of an IVOA service.
    """
    _test_get_access_url(service, capability, expected)

@unittest.mock.patch('requests.get', side_effect=_mocked_requests_get)
def _test_get_access_url(service, capability, expected, mock_get):
    orig_caps = commons.get_access_url.caps
    url = None
    try:
        # Force resource-caps to get called.
        commons.get_access_url.caps = {}
        url = commons.get_access_url(service, 'https://example.com/reg/resource-caps', capability)
        if url:
            assert url == expected
            tc = unittest.TestCase()
            expected_url = requests.utils.parse_url(expected)
            tc.assertIn(unittest.mock.call('https://example.com/reg/resource-caps'), mock_get.call_args_list, 'No resource-caps called.')
            tc.assertIn(unittest.mock.call(f'https://example.com{expected_url.path.replace("endpoint", "capabilities")}'), mock_get.call_args_list, 'No capabilities called.')
    except RuntimeError as runtime_error:
        # Should only happen if the provided service is empty.
        assert url is None
        assert str(runtime_error) == f'No or invalid service provided ({service}).'
    finally:
        commons.get_access_url.caps = orig_caps

