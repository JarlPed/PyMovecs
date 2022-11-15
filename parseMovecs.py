# A Python replacement for NwChem's mov2asc program needed to convert NwChem's
#
# The structure of .movecs file seems to consist of records, each being:
#  record = [4-bytes-integer-size] [data bytes] [4-bytes-integer-size]
#
# (c) Tymofii Nikolaienko, 2017
#
# This file is a part of the JANPA project. 
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 
# 3. All advertising and/or published materials mentioning features or use 
#    of this software must display the following acknowledgement:
# 
#       This product includes components from JANPA package of programs
#       ( http://janpa.sourceforge.net/ ) developed by Tymofii Nikolaienko
# 
# 4. Neither the name of the developer, Tymofii Nikolaienko,  nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
# 
# 5. In case if the code of JANPA package code and/or its parts and/or any data 
#    produced with JANPA package of programs are published, the following citation 
#    for the JANPA package of programs should be given:
# 
#     T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem. (2014)
#     V. 1050, P. 15-22, DOI: 10.1016/j.comptc.2014.10.002
# 
# THIS SOFTWARE IS PROVIDED BY ''AS IS'' AND ANY EXPRESS OR IMPLIED 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL Tymofii Nikolaienko BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

import sys
import struct

from datetime import datetime

import numpy as np
#from pyquante2.basis.basisset import basisset



        
class parseMovecs:
    
    
    
    def __init__(self, fname):
        
        self.f = open(fname, 'rb') 
        self.basissum, self.geomsum, self.bqsum = ('', '', '')
        self.extractFile()
        self.f.close()
        del(self.f)
        
    
    def read_rec(self, doPrint = False, convert = None):
        """ reads one record from .movecs file (global variable f) and converts 
            the read bytes if necessary """        
        #global f
        sz1 = struct.unpack("<L", self.f.read(4))[0]
        msg = self.f.read(sz1)
        sz2 = struct.unpack("<L", self.f.read(4))[0]
        if doPrint:
            print ('%d -> "%s" <- %d') % (sz1, msg, sz2)
            print ()
        if convert is None:
            return msg
        else:
            return struct.unpack(convert, msg)[0]
        
    def ret_double_arr(self):
        """ Reads a series of doubles, determining the number of values to
            read as data_size/8 """
        ret = []
        double_arr = self.read_rec()
        i = 0
        while i < len(double_arr):
            ret.append( struct.unpack("<d", double_arr[i:i+8])[0] )
            i += 8
        return ret
    
    def ret_int64_arr(self):
        """ Reads a series of doubles (that really are integers), determining the number of values to
            read as data_size/8 """
        ret = []
        double_arr = self.read_rec()
        i = 0
        while i < len(double_arr):
            ret.append( int(  struct.unpack("<d", double_arr[i:i+8])[0] )  )
            i += 8
        return ret
    
    def read_rec_nMO(self, doPrint = False, convert = None):
        """ reads one record from .movecs file (global variable f) and converts 
            the read bytes if necessary """        
        #global f
        sz1 = struct.unpack("<L", self.f.read(4))[0]
        msg = self.f.read(sz1)
        sz2 = struct.unpack("<L", self.f.read(4))[0]

        if convert is None:
            return msg
        else:
            return [ struct.unpack(convert, msg[j*8:j*8+8] )[0]  for j in range( int( len(msg)/8 ) )   ]
    
    
    
    
    def extractFile(self):
    
        hashsums = [self.basissum, self.geomsum, self.bqsum ]
    
        id = self.read_rec()
        i = 0
        for _ in range(3):      # basissum, geomsum, bqsum  // all in 32 plain characters
            hashsums[_] =   id[i:i+32].decode('utf-8')
            i += 32        
        self.me = id[i:i+20].decode('utf-8').strip()
    
    
    
        stringCtime = str(id[i+20:], "utf-8").strip()
        self.CreationTime = datetime.strptime(stringCtime, "%c")
        

    
        self.method = self.read_rec().decode('utf-8').strip()  # scftype20
    
        self.lentit =  self.read_rec(convert = "<Q")
    
        self.title = self.read_rec()  # title 
    
        self.lenbas =  self.read_rec(convert = "<Q") 
    
        sect = self.read_rec()
        self.basisName =  sect.decode('utf-8')   # basis_name
    
        self.nsets = self.read_rec(convert = "<Q")          # nsets
        self.nbf = self.read_rec(convert = "<Q")         # NBF Number of basis functions
        
    
            
        self.nMOs = self.read_rec_nMO(convert = "<Q")    # nmo(i) -- perhaps, the number if MOs in each set
    
        self.Sets = []
    
        for jset in range(self.nsets):
            temp = {}
            
            temp["OccupationNumbers"] =  self.ret_int64_arr()  
            
            temp["EigenValues"] =  self.ret_double_arr()   # Eigenvalues
            temp["EigenVectors"] = []
            for _ in range(self.nMOs[jset]): # Eigenvectors - the numerical value of the vave function phi. 
                temp["EigenVectors"].append(  self.ret_double_arr() )
            
            self.Sets.append(temp)
            
    
        # additional two doubles
        self.TotalEnergy, self.RepulsionEnergy =  self.ret_double_arr()     # energy, enrep
    

    
    def densityMatrix(self, setN=0):
        C = np.matrix(self.Sets[setN]["EigenVectors"] ).transpose() # shape = nMOs , nbf
        P = np.ones((self.nbf, self.nbf  ))
        for mu in range(self.nbf) :
            for nu in range(self.nbf ):
                #P[mu][nu] = 2 * sum(  [   C[mu,a]*np.conjugate(C[nu, a])   for a in range(self.nbf)   ])
                P[mu][nu] = 2 * sum(  [   C[mu,a]*C[nu, a]   for a in range(self.nMOs[0])   ])
        return P


if __name__ == "__main__":
    
    
    
    if len(sys.argv) < 2:
        sys.argv.append('./scfH2_sto-3g.movecs')
        
        #print ("Usage: python mov2asc.py filename.movecs" )
        #quit()

        
    fname = sys.argv[1]
    parser = parseMovecs(fname)
    print(parser.densityMatrix())

