# External imports
import os
import warnings
import requests
import urllib.parse
import tarfile
import shutil
import argparse

# Internal imports
from galapy.configuration import rcParams
from galapy.configuration.filesystem import _find_or_create_root_dir
import galapy.internal.globs as GP_GBL

################################################################################

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

################################################################################

def download_file ( url, out_file, overwrite = False, verbose = True, get_kw = {} ) :
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
        response = requests.get(url, **get_kw)
        response.raise_for_status()
    except requests.exceptions.HTTPError as http_err :
        print( f'HTTP error occurred: {http_err}' )
    except Exception as err :
        print( f'Other error occurred: {err}' )
        raise
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
            raise
        else :
            if verbose : print( f"Successfully created the directory {out_dir}" )

    return open(out_file, "wb").write(_bytes)

################################################################################

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
                    url = GP_GBL.DATA_URL.format( '/'.join( [ *path_line, file_name ] ) ),
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

################################################################################
# def search_file () :
#     pass
################################################################################

def download_database ( loc = None, name = None, version = None,
                        database = None, overwrite = False, verbose = True ) :
    """ A function for downloading the database from a github repository.
    (Or from any URL that can be converted into JSON and contains a 'tarball_url'
     attribute with the complete url to the tarball to be downloaded)

    Parameters:
    -----------
    loc : str
        Location in the filesystem where to download and extract the database.
        Default is None and will download in the default position specified in
        galapy.rcParams
    name : str
        Name to assign to the directory inside loc where the database will be extracted
        Default is None and will be set to the default name in galapy.internal.globs.DATA_DIR
    version : str
        The database version to download (Default is None and will download the latest version)
    database : str
        A valid URL to a database. Typically this would be in the form
        ``'https://api.github.com/repos/<user>/<repository>/releases'``.
        Any URL that can be converted by function requests.request.json() to a dictionary
        containing a valid URL to a tarball associated to key 'tarball_url' is valid.
    overwrite : bool
        Whether to overwrite whatever present at position ``os.path.join(loc, name)``.
        Default is ``False``, changing this SHOULD BE DONE AT OWN RISK.
    verbose : bool
        Whether to download verbosely. Default is ``True``

    Returns:
    --------
    : None
    """
    # a future implementation should give the possibility to add a new location to
    # the rcParams permanently
    from json import JSONDecodeError

    if loc is None :
        loc = rcParams['datapath'][-1]
    else :
        # checks the path already exists:
        if not os.path.exists( loc ) :
            # first, if the path is not present try to make a new directory:
            try :
                os.mkdir( loc )
            except OSError :
                if not os.path.isdir( loc ) :
                    raise
        elif not os.path.isdir( loc ) :
            # if the path exists but it is not a directory raise
            raise OSError( f'Intended galapy root directory {loc} is actually a file' )

        # return absolute path to selected location
        loc = os.path.abspath( loc )
        warnings.warn( f'The database will be downloaded in {loc}. '
                       'To add the custom location permanently to your '
                       'search paths modify the file '
                       f'{os.path.join(_find_or_create_root_dir(), "conf.ini")}' )
        
    if version is None :
        version = 'latest'
    elif version == 'latest' :
        pass
    elif version == 'default' :
        version = f'v{DATA_VERSION}'
    else :
        version = f'v{version}'

    if name is None :
        name = GP_GBL.DATA_DIR        
    destpath = os.path.join( loc, name )
    # If overwrite is set to False
    # and the path already exists leave the file as is
    # and return
    if os.path.exists( destpath ) and os.path.isdir( destpath ) :
        if not overwrite :
            warnings.warn( f"Directory {destpath} already exists." )
            return;
        else :
            shutil.rmtree( destpath )
    elif os.path.exists( destpath ) and not os.path.isdir( destpath ) :
        if not overwrite :
            warnings.warn( f"{destpath} is an already existing File." )
            return;
        else :
            os.remove( destpath )    

    # Check on input variable database
    if database is None :
        database = GP_GBL.DATABASE
    elif not is_url( database ) :
        raise OSError( f"The database provided {database:s} is not a valid URL." )

    #######################################################################
    # Download database

    # Request for a given release and check the URL is correct and working
    release = requests.get( database )
    try :
        release.raise_for_status()
    except requests.exceptions.HTTPError as http_err :
        print( f'HTTP error occurred: {http_err}' )
        raise
    except Exception as err :
        print( f'Other error occurred: {err}' )
        raise

    # converting json from list to dict
    releaselist = release.json()
    releasedict = {'latest' : releaselist[0]}
    releasedict.update(
        { rel['tag_name'] : rel for rel in releaselist }
    )

    # Extract from JSON the tarball URL and check for exceptions
    try :
        tarurl = releasedict[version]['tarball_url']
    except JSONDecodeError as err :
        print( f'Provided URL cannot be converted to JSON dictionary' )
        raise
    except Exception as err :
        print( f'An error occurred: {err}' )
        raise

    # Try downloading the tarball and eventually raise exception.
    where = os.path.join( loc, f'{name}.tar.gz' )
    try :
        if rcParams[ 'verbose_downloads' ] :
            print( 'downloading tarball, this might require some time ...' )
        nbytes = download_file( tarurl, where,
                                overwrite = overwrite,
                                verbose = verbose )
    except Exception as err :
        print( f'An error occurred: {err}' )
        raise
        
    if rcParams[ 'verbose_downloads' ] : print( f'Downloaded {nbytes} bytes' )

    #######################################################################
    # Extract and rename the archive
    
    # Open tarball
    archive = tarfile.open(where)

    # Get common prefix name (e.g. galapy_database-1.0.0)
    archive_root = os.path.commonprefix( archive.getnames() )

    # Extract everything from tarball on assigned position
    archive.extractall( loc )

    # Close tarball
    archive.close()

    # Rename extracted folder to assigned name
    os.rename(
        os.path.join( loc, archive_root ), # previous name
        destpath                           # new name
    )
    
    return;

################################################################################

def _entrypoint_download_database () :

    ####################################################################
    # Read command-line arguments:
    
    parser = argparse.ArgumentParser( description = 'options' )
    parser.add_argument( '--location', '-l',
                         dest = 'loc',
                         type = str,
                         default = None,
                         help = (
                             'Location in the filesystem where to download and ' +
                             'extract the database. ' +
                             f'DEFAULT: {rcParams["datapath"][-1]}'
                         ) )
    parser.add_argument( '--name', '-n',
                         dest = 'name',
                         type = str,
                         default = 'galapy_database',
                         help = (
                             'Name to assign to the directory inside ``loc`` where ' +
                             'the database will be extracted. ' +
                             'DEFAULT: "galapy_database"'
                         ) )
    parser.add_argument( '--version', '-v',
                         dest = 'version',
                         type = str,
                         default = 'latest',
                         help = (
                             'The database version to download. ' +
                             'DEFAULT: "latest"'
                         ) )
    args = parser.parse_args()
    
    ####################################################################

    return download_database(
        loc = args.loc,
        name = args.name,
        version = args.version,
        database = GP_GBL.DATABASE,
        overwrite = False,
        verbose = True
    )

################################################################################
