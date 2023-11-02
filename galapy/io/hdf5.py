import os, h5py, numpy, warnings
from collections.abc import MutableMapping as MM

def recursive_create ( group, indict ) :
    """Dumps an input dictionary to a HDF5 file
    """
    for k, v in indict.items() :
        try :
            if isinstance( v, MM ) :
                new_group = group.create_group( k )
                recursive_create( new_group, v )
            else :
                if hasattr( v, '__len__' ) and not isinstance( v, str ) :
                    group.create_dataset( k, data = v )
                else :
                    if v is not None :
                        group.attrs[k] = v
                    else :
                        group.attrs[k] = 'none'
        except Exception as err :
            print( f"- Skipping field '{k}' as it raises "
                   f"exception {type(err).__name__} with message:\n  '{err}'" )
    return;

def recursive_read ( group, outdict ) :
    """Reads a nested dictionary from a HDF5 file
    """
    outdict.update( { k : v if v != 'none' else None for k,v in dict(group.attrs).items()} )
    for k, v in group.items() :
        try :
            if isinstance( v, h5py.Dataset ) :
                outdict[k] = group.get(k)[:]
            if isinstance( v, h5py.Group ) :
                outdict[k] = {}
                recursive_read( v, outdict[k] )
        except Exception as err :
            print( f"- Skipping dataset '{k}' as it raises "
                   f"exception {type(err).__name__} with message:\n  '{err}'" )
            
    return;

def write_to_hdf5 ( outfile, metadata, hard = False, **groups ) :
    """Automatic dump of several dictionaries into hdf5 file as separate groups in file.
    Also support generic metadata
    
    Parameters
    ----------
    outfile : string
        /path/to/outputfile.hdf5 destination file. The path should exist on the filesystem
    metadata : dict
        dictionary with global metadata
    hard : bool
        (Optional, default=False) if True eventually overwrites an already existing
        file with the same name
    **groups : Keyword arguments
        (Optional) Dictionaries to save to the file. Each dictionary will be saved to the
        output file as a group with its corresponding keyword name as keyword.

    Returns
    -------
    None
    
    See Also
    --------
    load_from_hdf5 : loads all the groups in a HDF5 file
    """

    # check for overwrite
    if not hard and os.path.isfile( outfile ) :
        raise IOError( f"file\n{outfile}\nalready exist, "
                              "to overwrite it set argument hard=True." )

    # check directory exists
    outdir = os.path.dirname(outfile)
    if len( outdir ) == 0 :
        outdir = os.path.curdir
        outfile = os.path.join(outdir,outfile)
    elif not os.path.isdir( outdir ) :
        raise IOError( 'the provided path {os.path.dirname(outfile)} does not exist' )

    with h5py.File( outfile, "w" ) as f :
            
        # store meta-data
        for k,v in metadata.items() :
            f.attrs[k] = v

        # store groups:
        for name, dictionary in groups.items() :
            group = f.create_group(name)
            recursive_create( group, dictionary )
        
    print( f"Wrote on file {outfile}" )
    return;

def load_from_hdf5 ( infile ) :
    """Loads the groups in a hdf5 file into a dictionary.

    Parameters
    ----------
    infile : string
        /path/to/inputfile.hdf5 input file
    
    Returns
    -------
    : dict
        nested dictionary 
    
    See Also
    --------
    write_to_hdf5 : dump of several dictionaries into hdf5 file as separate groups in file
    
    Warning
    -------
    This is a simplified reader intended to work with hdf5 files generate
    using the ``galapy.io.hdf5.write_to_hdf5`` method.
    It is not guaranteed to work with any HDF5. 
    It assumes that at the top level all the fields can be cast to dictionaries.
    If this is not possible it raises a warning.
    """

    # check if file exists
    if not os.path.isfile( infile ) :
        raise IOError( f"file {infile} does not exist." )

    outdict = {}
    with h5py.File( infile, "r" ) as f :
        
        # load meta-data
        outdict['metadata'] = dict(f.attrs)

        # load groups
        for name, todict in f.items() :
            outdict[name] = {}
            recursive_read( todict, outdict[name] )
            
    print( f"Read from file {infile}" )
    return outdict
