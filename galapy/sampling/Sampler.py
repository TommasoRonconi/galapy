#############################################################################################

# External imports
import numpy
import sys

#############################################################################################

class Sampler () :
    
    def __init__ ( self, loglikelihood, ndim, 
                   sampler = 'dynesty',
                   prior_transform = None,
                   nwalkers = None,
                   pool = None,
                   dynesty_sampler_kw = {}, 
                   emcee_sampler_kw = {} ) :
        """
        """
        
        self.sampler = None
        self.which_sampler = sampler
        
        ######################################################################################
        # Dynamic Nested-Sampling with dynesty
        if self.which_sampler == 'dynesty' :
            if prior_transform is None or not hasattr(prior_transform, '__call__'):
                raise RuntimeError( 'You have to pass a function to `prior_transform` '
                                    'when running with the dynesty sampler.')
                
            from dynesty import DynamicNestedSampler
            self.sampler = DynamicNestedSampler( 
                loglikelihood = loglikelihood,
                prior_transform = prior_transform,
                ndim = ndim,
                pool = pool,
                **dynesty_sampler_kw
            )

        ######################################################################################
        # MCMC Ensemble sampling with emcee
        if self.which_sampler == 'emcee' :
            if nwalkers is None :
                raise RuntimeError( 'You have to pass an integer number of walkers to `nwalkers` '
                                    'when running with the emcee sampler.')
            from emcee import EnsembleSampler
            self.sampler = EnsembleSampler( nwalkers = nwalkers,
                                            ndim = ndim, 
                                            log_prob_fn = loglikelihood, 
                                            pool = pool, 
                                            **emcee_sampler_kw )
        
        ######################################################################################
        # Maybe in the future we'll add more sampling options ...    
        if self.sampler is None :
            raise AttributeError( f'The sampler chosen "{self.which_sampler}" is not valid.'
                                  'Valid samplers are ["dynesty", "emcee"].' )
            
    def run_sampling ( self, pos = None, nsample = None,
                       dynesty_sampling_kw = {}, emcee_sampling_kw = {} ) :
        """
        """
        from time import time
        
        if self.which_sampler == 'dynesty' :
            tstart = time()
            self.sampler.run_nested( **dynesty_sampling_kw )
            ndur = time() - tstart
            print( f'\nDone dynesty (dynamic) in {ndur} seconds' )
            # Print summary of the results:
            self.sampler.results.summary()
            
        if self.which_sampler == 'emcee' :
            if pos is None :
                raise RuntimeError( '`pos = None` but you should provide an initial position '
                                    'for the walkers when running with emcee.')
            if nsample is None :
                raise RuntimeError( '`nsample = None` but you should provide a number '
                                    'of samples to draw when running with emcee.')
            tstart = time()
            self.emcee_state = self.sampler.run_mcmc( pos, nsample, **emcee_sampling_kw )
            ndur = time() - tstart
            print( f'\nDone emcee in {ndur} seconds' )
            
        return;

    def save_results ( self, out_dir = '', name = '', pickle_sampler = False ) :
        """ Store the results of the last sampling run.
        NOTE THAT the type of output depends on the sampler used.
        
        Parameters
        ----------
        out_dir : string
          Position in the filesystem where the results will be stored.
          Default to the directory where the command has been called.
        name : string 
          A string identifying the run that will be saved. 
          By default it will use a string with the current date+time 
        pickle_sampler : bool
          If set to `True` the sampler will be pickled for future re-use.
        
        Returns
        -------
        : None
        """
        import os
        import pickle

        # If no output directory is passed, set it to the current working directory
        if len( out_dir ) == 0 or out_dir is None :
            out_dir = os.getcwd()

        # First check whether the required output directory exists,
        # if not it will be created in the correct position of the file-system
        if not os.path.isdir( out_dir ) :
            try :
                # creates multi-level subdirs. similarly to the *Nix command `mkdir -p`
                # (while os.mkdir() only allows to create the highest level directory)
                os.makedirs( out_dir ) 
            except OSError:
                print ( f"Creation of the directory {out_dir} failed" )

        # If no name for the current run has been provided, set it to
        # a string with the current date+time: 'year+month+day+hour+minute'
        if len(name) == 0 :
            from time import localtime as lt
            now = lt()
            name = ( f'{now.tm_year:04}{now.tm_mon:02}{now.tm_mday:02}' +
                     f'{now.tm_hour:02}{now.tm_min:02}' )

        # Set the string with the output base name
        outbase = os.path.join( out_dir, name )
            
        if self.which_sampler == 'dynesty' :
            with open( '_'.join( [ outbase, 'dynesty_results.pickle' ] ), 'wb' ) as f :
                pickle.dump( self.sampler.results, f )
            if pickle_sampler :
                with open( '_'.join( [ outbase, 'dynesty_sampler.pickle' ] ), 'wb' ) as f :
                    pickle.dump( self.sampler, f )

        if self.which_sampler == 'emcee' :
            # Pickling sampler-state:
            with open( '_'.join( [ outbase, 'emcee_state.pickle' ] ), 'wb' ) as f :
                pickle.dump( self.emcee_state, f )
            # Storing the flattened chains:
            numpy.save( '_'.join( [ outbase, 'emcee_flatchain' ] ),
                        self.sampler.flatchain )
            # Storing the flattened probabilities:
            numpy.save( '_'.join( [ outbase, 'emcee_flatlnprob' ] ),
                        self.sampler.flatlnprobability )
            # Storing the acceptance-fraction:
            numpy.save( '_'.join( [ outbase, 'emcee_axeptfrac' ] ),
                        self.sampler.acceptance_fraction )
            if pickle_sampler :
                with open( '_'.join( [ outbase, 'emcee_sampler.pickle' ] ), 'wb' ) as f :
                    pickle.dump( self.sampler, f )

        print( f'Results stored in files with prefix: {outbase}' )

        return;

#############################################################################################
