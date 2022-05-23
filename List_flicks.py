
"""
Ordered list of flicks files in a directory
"""

import os
flicks_files = [f.replace('flicks.','') for f in os.listdir("./") if os.path.isfile(os.path.join("./", f)) and 'flicks.' in f and len(f.split('.'))==2 and f.split('.')[1].isnumeric()]
flicks_files.sort()
print(flicks_files)
