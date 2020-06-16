#!C:\Users\sth19\PycharmProjects\Genome_Graph_project\pyGFA\venv\Scripts\python.exe
# EASY-INSTALL-ENTRY-SCRIPT: 'coverage==5.1','console_scripts','coverage3'
__requires__ = 'coverage==5.1'
import re
import sys
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(
        load_entry_point('coverage==5.1', 'console_scripts', 'coverage3')()
    )
