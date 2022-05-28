# External includes
import os

def _find_home():
    """Locates and return the home directory (or best approximation) on this
    system.

    Raises
    ------
    OSError
        If the home directory cannot be located - usually means you are running
        GalaPy on some obscure platform that doesn't have standard home
        directories.

    Notes
    -----
    This routine has been adapted from the Astropy Collaboration (http://www.astropy.org)
    """
    try:
        homedir = os.path.expanduser('~')
    except Exception:
        # Linux, Unix, AIX, OS X
        if os.name == 'posix':
            if 'HOME' in os.environ:
                homedir = os.environ['HOME']
            else:
                raise OSError('Could not find unix home directory to search for '
                              'galapy config dir')
        elif os.name == 'nt':  # This is for all modern Windows (NT or after)
            if 'MSYSTEM' in os.environ and os.environ.get('HOME'):
                # Likely using an msys shell; use whatever it is using for its
                # $HOME directory
                homedir = os.environ['HOME']
            # See if there's a local home
            elif 'HOMEDRIVE' in os.environ and 'HOMEPATH' in os.environ:
                homedir = os.path.join(os.environ['HOMEDRIVE'],
                                       os.environ['HOMEPATH'])
            # Maybe a user profile?
            elif 'USERPROFILE' in os.environ:
                homedir = os.path.join(os.environ['USERPROFILE'])
            else:
                try:
                    import winreg as wreg
                    shell_folders = r'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders'  # noqa: E501
                    key = wreg.OpenKey(wreg.HKEY_CURRENT_USER, shell_folders)

                    homedir = wreg.QueryValueEx(key, 'Personal')[0]
                    key.Close()
                except Exception:
                    # As a final possible resort, see if HOME is present
                    if 'HOME' in os.environ:
                        homedir = os.environ['HOME']
                    else:
                        raise OSError('Could not find windows home directory to '
                                      'search for galapy config dir')
        else:
            # for other platforms, try HOME, although it probably isn't there
            if 'HOME' in os.environ:
                homedir = os.environ['HOME']
            else:
                raise OSError('Could not find a home directory to search for '
                              'galapy config dir - are you on an unsupported '
                              'platform?')
    return homedir

def _find_or_create_root_dir() :
    """ Search for the GalaPy rootdir in user's home directory.
    If it cannot find it, one is created in ``<home>/.galapy``.
    """

    rootdir = os.path.join( _find_home(), '.galapy' )

    # checks the path already exists:
    if not os.path.exists( rootdir ) :
        # first, if the path is not present try to make a new directory:
        try :
            os.mkdir( rootdir )
        except OSError :
            if not os.path.isdir( rootdir ) :
                raise
    elif not os.path.isdir( rootdir ) :
        # if the path exists but it is not a directory raise
        raise OSError( f'Intended galapy root directory {rootdir} is actually a file' )

    # return absolute path to root directory
    return os.path.abspath( rootdir )
