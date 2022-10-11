#############################################################################################
# External imports

import numpy
import sys

#############################################################################################
# Internal imports

from galapy.internal.utils import now_string

#############################################################################################

emcee_default_sampler_kw = { 'moves' : None, 'args' : None, 'kwargs' : None }
emcee_default_sampling_kw = { 'log_prob0' : None, 'rstate0' : None, 'blobs0' : None,
                              'tune' : False, 'skip_initial_state_check' : False,
                              'thin_by' : 1, 'thin' : None, 'store' : True,
                              'progress' : True, 'progress_kwargs' : None }

#############################################################################################

dynesty_default_sampler_kw = { 'bound' : 'multi', 'sample' : 'rwalk',
                               'update_interval' : 0.6, 'walks' : 25,
                               'bootstrap' : 0 }
dynesty_default_sampling_kw = { 'nlive_init' : 1024, 'maxiter_init' : None,
                                'maxcall_init' : None, 'dlogz_init' : 0.02,
                                'nlive_batch' : 1024, 'maxiter_batch' : None,
                                'maxcall_batch' : None,
                                'maxiter' : sys.maxsize, 'maxcall' : sys.maxsize,
                                'maxbatch' : sys.maxsize,
                                'n_effective' : None, 'print_progress' : True,
                                'stop_kwargs' : { 'target_n_effective' : 10000 } }

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

            kw = dynesty_default_sampler_kw
            kw.update( dynesty_sampler_kw )
            from dynesty import DynamicNestedSampler
            self.sampler = DynamicNestedSampler( 
                loglikelihood = loglikelihood,
                prior_transform = prior_transform,
                ndim = ndim,
                pool = pool,
                **kw
            )
            
        ######################################################################################
        # MCMC Ensemble sampling with emcee
        if self.which_sampler == 'emcee' :
            if nwalkers is None :
                raise RuntimeError( 'You have to pass an integer number of walkers to `nwalkers` '
                                    'when running with the emcee sampler.')
            kw = emcee_default_sampler_kw
            kw.update( emcee_sampler_kw )
            from emcee import EnsembleSampler
            self.sampler = EnsembleSampler( nwalkers = nwalkers,
                                            ndim = ndim, 
                                            log_prob_fn = loglikelihood, 
                                            pool = pool, 
                                            **kw )
        
        ######################################################################################
        # Maybe in the future we'll add more sampling options ...    
        if self.sampler is None :
            raise AttributeError( f'The sampler chosen "{self.which_sampler}" is not valid.'
                                  'Valid samplers are ["dynesty", "emcee"].' )
            
    def run_sampling ( self, pos = None, nsample = None,
                       dynesty_sampling_kw = {}, emcee_sampling_kw = {} ) :
        """ Run the sampling with the chosen sampler.
        
        Parameters
        ----------
        pos : iterable
        nsample : int
        
        Returns
        -------
        : None
        """
        from time import time
        
        if self.which_sampler == 'dynesty' :
            kw = dynesty_default_sampling_kw
            kw.update( dynesty_sampling_kw )
            tstart = time()
            self.sampler.run_nested( **kw )
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
            kw = emcee_default_sampling_kw
            kw.update( emcee_sampling_kw )
            tstart = time()
            self.emcee_state = self.sampler.run_mcmc( pos, nsample, **kw )
            ndur = time() - tstart
            print( f'\nDone emcee in {ndur} seconds' )
            
        return;

    def return_samples_logl_weights ( self ) :

        outtuple = None
        if self.which_sampler == 'dynesty' :
            weights = numpy.exp( self.sampler.results.logwt -
                                 self.sampler.results.logz[-1] )
            outtuple = ( self.sampler.results.samples,
                         self.sampler.results.logl,
                         weights )
        if self.which_sampler == 'emcee' :
            outtuple = ( self.sampler.flatchain,
                         self.sampler.flatlnprobability,
                         numpy.ones_like( self.sampler.flatlnprobability ) )
        if outtuple is None :
            raise RuntimeError(
                "Function must be called after sampling run"
            )
        
        return outtuple

    def save_results ( self, outbase = '',
                       pickle_sampler = False,
                       pickle_raw = True ) :
        """ Store the results of the last sampling run.
        NOTE THAT the type of output depends on the sampler used.
        
        Parameters
        ----------
        outbase : string
          Position in the filesystem where the results will be stored.
          Default to the directory where the command has been called with
          a current date+time as name.
        pickle_sampler : bool
          If set to `True` the sampler will be pickled for future re-use.
        pickle_raw : bool
          If set to `True` the raw results of the sampling run will be pickled.
        
        Returns
        -------
        : None
        """
        import os
        import pickle

                    
        # If no name for the current run has been provided, set it to
        # a string with the current date+time: 'year+month+day+hour+minute'
        if len(outbase) == 0 :
            outbase = now_string()
            
        if self.which_sampler == 'dynesty' :
            if pickle_raw :
                with open( '_'.join( [ outbase, 'dynesty_results.pickle' ] ), 'wb' ) as f :
                    pickle.dump( self.sampler.results, f )
            if pickle_sampler :
                with open( '_'.join( [ outbase, 'dynesty_sampler.pickle' ] ), 'wb' ) as f :
                    pickle.dump( self.sampler, f )

        if self.which_sampler == 'emcee' :
            if pickle_raw :
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

        print( f'Sampler results stored in files with prefix: {outbase}' )

        return;

#############################################################################################
