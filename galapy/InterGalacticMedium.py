# External imports
import numpy

# Internal imports
import galapy.internal.globs as GP_GBL
from galapy.internal.data import DataFile

class IGM () :
    
    LymanLimit = 912.

    def __init__ ( self ) :
        import os
        
        # Shorter version for convenience
        self.lLL = self.LymanLimit
        
        # Get parameters from file
        intuple = numpy.genfromtxt( DataFile( *GP_GBL.IGM_Inoue14 ).get_file(), 
                                    unpack = True )
        ( self.line, self.lambda_j, 
          self.ALAFj1, self.ALAFj2, self.ALAFj3, 
          self.ADLAj1, self.ADLAj2 ) = intuple
        self.j_lambda = 1. / self.lambda_j
        
    def LS_LAF ( self, lobs, zobs ) :
        
        opt_depth = numpy.zeros_like(lobs)
        
        for jj, (lj, jl) in enumerate(zip(self.lambda_j, 
                                          self.j_lambda)) :
            
            mask = (lj < lobs) & (lobs < lj*(1+zobs))
            l1 = lobs < 2.2*lj
            l2 = (2.2*lj <= lobs) & (lobs < 5.7*lj)
            l3 = 5.7*lj <= lobs
            m1,m2,m3 = mask&(l1,l2,l3)
            
            opt_depth[m1] += self.ALAFj1[jj]*(jl * lobs[m1])**1.2
            opt_depth[m2] += self.ALAFj2[jj]*(jl * lobs[m2])**3.7
            opt_depth[m3] += self.ALAFj3[jj]*(jl * lobs[m3])**5.5
            
        return opt_depth
    
    def LS_DLA ( self, lobs, zobs ) :
        
        opt_depth = numpy.zeros_like(lobs)
        
        for jj, (lj, jl) in enumerate(zip(self.lambda_j, 
                                          self.j_lambda)) :
            
            mask = (lj < lobs) & (lobs < lj*(1+zobs))
            l1 = lobs < 3.0*lj
            l2 = 3.0*lj <= lobs
            m1,m2 = mask&(l1,l2)
            
            opt_depth[m1] += self.ADLAj1[jj]*(jl * lobs[m1])**2
            opt_depth[m2] += self.ADLAj2[jj]*(jl * lobs[m2])**3
            
        return opt_depth

    def LC_LAF ( self, lobs, zobs ) :
    
        # preparing quantities
        opt_depth = numpy.zeros_like(lobs)
        mask = self.lLL < lobs
        lu = lobs/self.lLL

        # different behaviour depending on source redshift
        if zobs < 1.2 :

            # 1
            m1 = mask & (lobs < self.lLL*(1+zobs))
            opt_depth[m1] += 0.325 * (
                lu[m1]**1.2 - (1+zobs)**(-0.9)*lu[m1]**2.1
            )

            # 2 not necessary, return
            return opt_depth

        if (1.2 <= zobs)&(zobs < 4.7) :

            # 1
            m1 = mask & (lobs < 2.2*self.lLL)
            opt_depth[m1] += (
                2.55e-2 * (1+zobs)**1.6 * lu[m1]**2.1 +
                0.325 * lu[m1]**1.2 - 0.250 * lu[m1]**2.1
            )

            # 2
            m2 = mask & ((2.2*self.lLL <= lobs) & 
                         (lobs < self.lLL*(1+zobs)))
            opt_depth[m2] += 2.55e-2 * (
                (1+zobs)**1.6 * lu[m2]**2.1 - lu[m2]**3.7
            )

            # 3 not necessary, return
            return opt_depth

        if 4.7 <= zobs :

            # 1
            m1 = mask & (lobs < 2.2*self.lLL)
            opt_depth[m1] += (
                5.22e-4 * (1+zobs)**3.4 * lu[m1]**2.1 +
                0.325 * lu[m1]**1.2 - 3.14e-2 * lu[m1]**2.1
            )

            # 2
            m2 = mask & ((2.2*self.lLL <= lobs) & 
                         (lobs < 5.7*self.lLL))
            opt_depth[m2] += (
                5.22e-4 * (1+zobs)**3.4 * lu[m2]**2.1 +
                0.218 * lu[m2]**2.1 - 2.55e-2 * lu[m2]**3.7
            )

            # 3
            m3 = mask & ((5.7*self.lLL <= lobs) & 
                         (lobs < self.lLL*(1+zobs)))
            opt_depth[m3] += 5.22e-4 * (
                (1+zobs)**3.4 * lu[m3]**2.1 - lu[m3]**5.5
            )

            # 4 not necessary, return
            return opt_depth

        # this should never be raised
        raise Exception
        
    def LC_DLA ( self, lobs, zobs ) :
    
        # preparing quantities
        opt_depth = numpy.zeros_like(lobs)
        mask = self.lLL < lobs
        lu = lobs/self.lLL

        # different behaviour depending on source redshift
        if zobs < 2.2 :

            # 1
            m1 = mask & (lobs < self.lLL*(1+zobs))
            opt_depth[m1] += (
                0.211 * (1+zobs)**2 - 
                7.66e-2 * (1+zobs)**2.3 * lu[m1]**(-0.3) -
                0.135 * lu[m1]**2
            )

            # 2 not necessary, return
            return opt_depth

        else :

            # 1
            m1 = mask & (lobs < 3.0*self.lLL)
            opt_depth[m1] += (
                0.634 + 4.7e-2 * (1+zobs)**3 -
                1.78e-2 * (1+zobs)**3.3 * lu[m1]**(-0.3) -
                0.135 * lu[m1]**2 - 0.291 * lu[m1]**(-0.3)
            )

            # 2
            m2 = mask & ((3.0*self.lLL <= lobs) & 
                         (lobs < self.lLL*(1+zobs)))
            opt_depth[m2] += (
                4.7e-2 * (1+zobs)**3 - 
                1.78e-2 * (1+zobs)**3.3 * lu[m2]**(-0.3) -
                2.92e-2 * lu[m2]**3
            )

            # 3 not necessary, return
            return opt_depth

        # this should never be raised
        raise Exception

    def optical_depth ( self, lobs, zobs ) :
        return (
            self.LS_LAF( lobs, zobs ) +
            self.LC_LAF( lobs, zobs ) + 
            self.LS_DLA( lobs, zobs ) +
            self.LC_DLA( lobs, zobs )
        )
    
    def transmission ( self, lobs, zobs ) :
        return numpy.exp( -self.optical_depth( lobs, zobs ) )
