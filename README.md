## ASOT, a suite of analysis routines for the ARMS code

> *"Or to take arms against **A** **S**ea **O**f **T**roubles"* 

    -- W. Shakespeare

  

Most of the routines are written in Python3. **To make them work, you must change the path near the top of each file, which by default is `'/Change/This/Path'` to the directory where the file `ASOT_Functions_Python.py` is kept.**

The structure of most Python files is as follows:<br>
    `frequently changed variables`<br><br>

    `import sys`<br>
    `sys.path[:0]=['/Change/This/Path']`<br>
    `from ARMS_ASOT_Functions import *`<br>
    `import matplotlib.pyplot as plt`<br><br>

    `routine specific code`
<br><br><br><br>
There are also routines to analyze [**QSL Squasher**](https://arxiv.org/abs/1609.00724), which is [hosted here](https://bitbucket.org/tassev/qsl_squasher).<br>

Please acknowledge the use of this software when publishing work which made significant use of it.<br>


Note that this repository is an effort to share a personal set of routines which have many embarrassingly poor optimizations and the like. Any input is appreciated.
