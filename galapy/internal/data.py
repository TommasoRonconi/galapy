# External imports
import os
import warnings
import requests
import urllib.parse

# Internal imports
from galapy.configuration import rcParams
from galapy.internal.globs import *

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

    Notes
    -----
    This routine has been adapted from the Astropy Collaboration (http://www.astropy.org)
    """
    u = urllib.parse.urlparse(string)
    # we can't just check that url.scheme is not an empty string, because
    # file paths in windows would return a non-empty scheme (e.g. e:\\
    # returns 'e').
    return u.scheme.lower() in ['http', 'https', 'ftp', 'sftp', 'ssh', 'file']

def download_file ( url, out_file, overwrite = False, verbose = True ) :
    """ Download file from known url into defined position

    Parameters
    ----------
    url : str
      The complete url of the file to download
    out_file : str
      Absolute path to output file
    overwrite : bool
      If set to True, the eventually already existing out_file is overwritten
    verbose : bool
      If set to True, prints additional informations while processing

    Returns
    -------
    : int
    size in bytes written to the out_file

    Raises
    ------
    OSError 
      If the input url string is not a valid URL. 
      If the output file directory does not exist and the system is not able to create it.
    requests.exceptions.HTTPError
      If an HTTP error occurred.
    Exception
      If some other error occurs.
    """

    if not is_url( url ) :
        raise OSError( f"The string provided {url:s} is not a valid URL." )

    # Extracting bytes from the provided URL:
    _bytes = b''
    try :
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.HTTPError as http_err :
        print( f'HTTP error occurred: {http_err}' )
    except Exception as err :
        print( f'Other error occurred: {err}' )
    else :
        _bytes = response.content

    # If overwrite is set to False
    # and the path already exist leave the file as is
    # returns zero (number of bytes written)
    if os.path.exists( out_file ) and not overwrite :
        warnings.warn( f"File {out_file} already exists." )
        return 0
        
    # Checking whether the directory for file output is accessible 
    out_dir = os.path.dirname( out_file )
    if not os.path.isdir( out_dir ) :
        if verbose : print( f"Directory {out_dir} does not exist, it will be created." )
        try :
            # creates multi-level subdirs. similarly to the *Nix command `mkdir -p`
            # (while os.mkdir() only allows to create the highest level directory)
            os.makedirs( out_dir ) 
        except OSError:
            print ( f"Creation of the directory {out_dir} failed" )
        else :
            if verbose : print( f"Successfully created the directory {out_dir}" )

    return open(out_file, "wb").write(_bytes)

class DataFile () :
    """
    - search in the rcParams[ 'datapath' ] list
    - if not present, download and store it 
    """

    def __init__ ( self, file_name, path_line ) :

        # allocate variable to store the file-path
        self._filepath = None

        # search in the available data-paths for the requested file
        for datapath in rcParams[ 'datapath' ] :
            _testpath = os.path.join( datapath, *path_line, file_name )
            if not os.path.exists( _testpath ) :
                continue
            elif os.path.isdir( _testpath ) :
                raise OSError( f'The DataFile {file_name} '
                               f'in position {os.path.join(datapath, *path_line)} '
                               'is a directory.' )
            else :
                self._filepath = _testpath

        # if the data-file is not present in the system download it from remote
        if self._filepath is None :
            where = os.path.join( rcParams[ 'datapath' ][ 0 ],
                                  *path_line,
                                  file_name )
            try :
                byte_count = download_file(
                    url = DATA_URL.format( '/'.join( [ *path_line, file_name ] ) ),
                    out_file = where,
                    overwrite = rcParams[ 'force_downloads' ],
                    verbose = rcParams[ 'verbose_downloads' ]
                )
            except Exception as err :
                print( f'An error occurred: {err}' )
                raise
            else :
                if rcParams[ 'verbose_downloads' ] : print( f'Downloaded {byte_count} bytes.' )
                self._filepath = where
            

    def get_file ( self ) :
        """ 
        - return the string with the cross-platform path
        """
        return self._filepath

# def search_file () :
#     pass

# def download_database ( loc = None ) :
#     pass
