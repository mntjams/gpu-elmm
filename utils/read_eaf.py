#!/usr/bin/python3
import numpy as np

scalar_flag = 1001

vector_flag = 1002
   

class eaf_reader:
  
    def __init__(self, base_name):
    
        self.base_name = base_name

        self.efh_name = base_name + ".efh"

        efh = open(self.efh_name, "rb")

        one = np.fromfile(efh, dtype="int32", count=1)[0]

        if (one != 1):
          print("Error, the data has the the opposite endianness!")
          raise ValueError
            
        
        float_size = np.fromfile(efh, dtype="int32", count=1)[0]
        
        if (float_size == 4):
            self.dtype = "float32"
        elif (float_size == 8):
            self.dtype = "float64"
        else:
            print("Unknown float size in " + efh_name)
            raise ValueError
          
        self.nxyz = np.fromfile(efh, dtype="int32", count=3).tolist()
        
        self.xs = np.fromfile(efh, dtype = self.dtype, count = self.nxyz[0])
        self.ys = np.fromfile(efh, dtype = self.dtype, count = self.nxyz[1])
        self.zs = np.fromfile(efh, dtype = self.dtype, count = self.nxyz[2])
        
        self.arrays = []
        #Now for simplicitly just assume that we know that the file contains velocity vectors only.
        #Therefore we skip reading the data description.
        self.arrays = self.arrays + [{"name":"u", "vector":True}]
        
    def read_frame(self, n_frame):
        eaf = open(self.base_name + "-" + str(n_frame) + ".eaf", "rb")
        
        data = []
        
        for array in self.arrays:
            if array["vector"]:
              array_data = np.fromfile(eaf, dtype=self.dtype).reshape([3]+self.nxyz, order='F')
            else:
              array_data = np.fromfile(eaf, dtype=self.dtype).reshape(self.nxyz, order='F')
              
            data = data + [{"name" : array["name"], "data" : array_data}]
            
        return data
   
    #returns the x-xoordinate of the ith index
    def x(self,i):
      return self.xs[i]
      
    #returns the y-xoordinate of the ith index
    def y(self,i):
      return self.ys[i]
      
    #returns the z-xoordinate of the ith index
    def z(self,i):
      return self.zs[i]
      
                           
if __name__ == '__main__':
   name = "frame-test"
   
   reader = eaf_reader(name)
   
   print(reader.x(0))
   print(reader.y(0))
   print(reader.z(0))
   
   data = reader.read_frame(0)
   
   print(data[0]["name"])
   
   u = data[0]["data"]
   
   print(np.shape(u))
   
   #first index for vectors: 0..x component, 1..y component, 2..z component
   
   #then x index, y index, z index
   print(u[0,0,0,0])
   print(u[1,1,0,0])
   print(u[2,0,0,1])
   
   
   
   
