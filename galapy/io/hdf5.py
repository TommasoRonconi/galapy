import os, h5py, numpy, warnings

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
    load_from_hdf5 : loads all the groups in a ``.hdf5`` file
    """

    # check right extension
    if outfile.split('.')[-1] != 'hdf5' :
        raise AttributeError( 'the provided file does not have a hdf5 extension' )

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
            for key, data in dictionary.items() :
                try :
                    group.create_dataset( key, data = data )
                except Exception as err :
                    print( f"- Skipping dataset '{key}' of group '{name}' as it raises "
                           f"exception {type(err).__name__} with message:\n  '{err}'" )
        
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
    using the ``scampy.io.hdf5.write_to_hdf5`` method.
    It is not guaranteed to work with any ``.hdf5``. 
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
            try :
                outdict[name] = { k : todict.get(k)[:] for k in todict.keys() }
            except Exception as err :
                warnings.warn(
                    f"- Skipping group '{name}' as it raises "
                    f"exception {type(err).__name__} with message:\n  '{err}'"
                )
            
    print( f"Read from file {infile}" )
    return outdict
