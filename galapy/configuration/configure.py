# External includes:
import configparser
import os

# Internal includes
from galapy.configuration.filesystem import _find_or_create_root_dir

_confstr = f"""
#############################
# GalaPy configuration file #
#############################

# Contains the default and user-defined paths
# to the input data-files and template used in the library
# Note that the system will search data-files in all the 
# paths listed below before trying to download the relevant
# files. Research will be made in order of appearence in this section.
[DATAPATHS]

DEFAULT={{0:s}}

# Add eventual new positions in the filesystem here:
# (the name given to the new position is not relevant)
# NEWPATH={os.path.join('absolute','path','to','new','position')}

# Here define boolean values that affect the generic behaviour 
# of the package. 
# It accepts values as true/false, yes/no, 1/0, on/off
[BOOLEANS]

# This variable decides whether to force overwriting of whatever present in a
# path that has been found to already exist during a download request. 
# Setting this to `true` is done at your own risk. 
FORCE_DOWNLOADS=false

# When set to true it will provide additional informations
# while downloading data-files
VERBOSE_DOWNLOADS=true
"""

def _find_or_create_config_file () :
    """ Searches the conf.ini file, it generates one if none exist

    Returns
    -------
    : str
      the absolute path to the configuration file
    
    Raises
    ------
    OSError
      when the configuration file is not present and it is
      not possible to generate any 
    """
    
    rootdir = _find_or_create_root_dir()
    conf_file = os.path.join( rootdir, 'conf.ini' )
    if not os.path.exists( conf_file ) :
        open( conf_file, "w" ).write( _confstr.format(rootdir) )
    elif os.path.isdir( conf_file ) :
        raise OSError( f"Intended galapy configuration file {conf_file}"
                       " is actually a directory" )

    return conf_file

def configure () :
    """ Finds/Generates the conf.ini file, then it reads it
    and converts it to a dictionary, effectively generating the
    rc-parameters.
    
    Returns
    -------
    : dict
      dictionary containing the configuration of the package
    """
    
    conf = configparser.ConfigParser()
    conf.read( _find_or_create_config_file() )

    ###########################################################
    # Use ConfigParser to generate running-command dictionary #
    ###########################################################

    rcdict = {}
    rcdict[ 'datapath' ] = [ v for _, v in conf.items( 'DATAPATHS' ) ]
    for k in conf[ 'BOOLEANS' ].keys() :
        rcdict[ k.lower() ] = conf[ 'BOOLEANS' ].getboolean( k )

    return rcdict
