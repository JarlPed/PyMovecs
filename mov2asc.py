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


#from pyquante2.basis.basisset import basisset



nbf = 1E5 # dummy value :  NBF Number of basis functions

        
        
def MOV2ASCII (fname):
    f = open(fname, 'rb') 
    
    def read_rec(doPrint = False, convert = None):
        """ reads one record from .movecs file (global variable f) and converts 
            the read bytes if necessary """        
        #global f
        sz1 = struct.unpack("<L", f.read(4))[0]
        msg = f.read(sz1)
        sz2 = struct.unpack("<L", f.read(4))[0]
        if doPrint:
            print ('%d -> "%s" <- %d') % (sz1, msg, sz2)
            print ()
        if convert is None:
            return msg
        else:
            return struct.unpack(convert, msg)[0]
        
    def print_double_arr():
        """ Reads a series of doubles, determining the number of values to
            read as data_size/8 """
        double_arr = read_rec()
        i = 0
        while i < len(double_arr):
            print ("%25.15E" % ( struct.unpack("<d", double_arr[i:i+8])[0] ) )
            i += 8
            if i % ((nbf*8)) == 0:
                print ()
        if (i % (nbf*8)) != 0:
            print ()
    
    def print_int64_arr():
        """ Reads a series of doubles (that really are integers), determining the number of values to
            read as data_size/8 """
        double_arr = read_rec()
        i = 0
        while i < len(double_arr):
            print ("%i" % ( struct.unpack("<d", double_arr[i:i+8])[0] ) )
            i += 8
            if i % ((nbf*8)) == 0:
                print ()
        if (i % (nbf*8)) != 0:
            print ()
    
    
    def read_rec_nMO(doPrint = False, convert = None):
        """ reads one record from .movecs file (global variable f) and converts 
            the read bytes if necessary """        
        #global f
        sz1 = struct.unpack("<L", f.read(4))[0]
        msg = f.read(sz1)
        sz2 = struct.unpack("<L", f.read(4))[0]
        if doPrint:
            print ('%d -> "%s" <- %d') % (sz1, msg, sz2)
            print ()
        if convert is None:
            return msg
        else:
            return [ struct.unpack(convert, msg[j*8:j*8+8] )[0]  for j in range( int( len(msg)/8 ) )   ]
    
    
    
    
    
    
    print ("# This is an NWChem movecs file translated by mov2asc")
    
    hashsums = ['basissum', 'geomsum', 'bqsum' ]
    
    id = read_rec()
    i = 0
    for _ in range(3):      # basissum, geomsum, bqsum  // all in 32 plain characters
        print ( hashsums[_] + ":\t" +  id[i:i+32].decode('utf-8')   )  
        i += 32        
    print(id[i:i+20].decode('utf-8').strip()) # scftype20
    
    
    
    stringCtime = str(id[i+20:], "utf-8").strip()
    CreationTime = datetime.strptime(stringCtime, "%c")
    print(CreationTime.ctime())  # date
    
    #print(id[i+20:].decode('utf-8'))  # date
    

    
    method = read_rec().decode('utf-8').strip()  # scftype20
    print(method)
    
    print ("lentit:\t%10d" % read_rec(convert = "<Q"))  # lentit
    
    title = read_rec()  # title 
    print ("Title:\t%s" % title.decode('utf-8'))
    
    print ("lenbas:\t%10d" % read_rec(convert = "<Q"))  # lenbas
    
    sect = read_rec()
    print (sect.decode('utf-8'))   # basis_name
    
    nsets = read_rec(convert = "<Q")          # nsets
    print ("nsets:\t %10d" % nsets  ) 
    nbf = read_rec(convert = "<Q")         # NBF Number of basis functions
    print ("Number of basis functions(NBF):\t%10d" % nbf)
    nMOs = []
    
            
    nMOs = read_rec_nMO(convert = "<Q")    # nmo(i) -- perhaps, the number if MOs in each set
    for j in range(len(  nMOs)):
        print ( "Number of MO's in set %i: %10d" %(j, nMOs[j] ) )
    
    for jset in range(nsets):
        print("Occupation numbers:")
        print_int64_arr()   # Occupation numbers - number of electrons occupying the orbital; for dft, nsets = 2 and Occupation numbers are  1; BUT there is actually 2 electrons occupying the same orbital. for regular scf; there is only a single set with occupation numbers always equal to 2.
        print("Eigenvalues:")
        print_double_arr()   # Eigenvalues
        print("Eigenvectors:")
        for _ in range(nMOs[jset]): # Eigenvectors - the numerical value of the vave function phi. 
            print_double_arr()
    
    # additional two doubles
    print("Total energy & Repulsive (nucleus-nucleus) energy")
    print_double_arr()     # energy, enrep
    
    
    #print '-'*80
    #sz = struct.unpack("<L", f.read(4))[0]
    #print 'next: ',sz
    

    
    


if __name__ == "__main__":
    
    
    
    if len(sys.argv) < 2:
        sys.argv.append('./dftH2_STO3G.movecs')
        
        #print ("Usage: python mov2asc.py filename.movecs" )
        #quit()

        
    fname = sys.argv[1]
    MOV2ASCII(fname)

