import os
import requests
import urllib.parse

# base_url = "https://raw.githubusercontent.com/TommasoRonconi/galapy/tree/main/galapy/internal/data"

def is_url ( string ) :
    """
    Test whether a string is a valid URL for :func:`download_file`.

    Parameters
    ----------
    string : str
        The string to test.

    Returns
    -------
    status : bool
        String is URL or not.

    Warning
    -------
    This function has been stolen from the Astropy Collaboration
    """
    u = urllib.parse.urlparse(string)
    # we can't just check that url.scheme is not an empty string, because
    # file paths in windows would return a non-empty scheme (e.g. e:\\
    # returns 'e').
    return u.scheme.lower() in ['http', 'https', 'ftp', 'sftp', 'ssh', 'file']

def download_file ( url, out_file, force = False ) :
    """ Download file from known url into defined position

    Parameters
    ----------
    url : str
      The complete url of the file to download
    out_file : str
      Absolute path to output file

    Returns
    -------
    : int
    size in bytes of the out_file
    """

    if not is_url( url ) :
        raise OSError( f"URL {url:s} does not link to any valid file." )
    
    r = requests.get( url )
    return open(out_file, "wb").write(r.content)

