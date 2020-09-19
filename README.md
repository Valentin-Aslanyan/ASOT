## ASOT, a suite of analysis routines for the ARMS code

> *"Or to take arms against **A** **S**ea **O**f **T**roubles"* 

    -- W. Shakespeare

  

Most of the routines are written in Python3. **To make them work, you must change the path near the top of each file, which by default is `'/Change/This/Path'` to the directory where the file `ASOT_Functions_Python.py` is kept.**

The structure of most Python files is as follows:
`frequently changed variables`
  
`import sys`
`sys.path[:0]=['/Change/This/Path']`
`from ARMS_ASOT_Functions import *`
`import matplotlib.pyplot as plt`
  
`routine specific code`
  
  
   
  
  
  

Note that this repository is an effort to share a personal set of routines which have many embarrassingly poor optimizations and the like. Any input is appreciated.
