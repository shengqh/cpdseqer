# In order to run the code directly for debug, I will need to run python cpdseqer/__main__py directly.
# Unfortunately, it will cause error:
#   File "cpdseqer/__main__.py", line 7, in <module>
#     from .__version__ import __version__
#   ModuleNotFoundError: No module named '__main__.__version__'; '__main__' is not a package
# We have to call the __main__.main from outside of package.
# So I create this debug.py to call the function

import re
import sys

from cpdseqer.__main__ import main

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
