import os
import webbrowser
from HTML import table


path = os.path.abspath(os.path.join("../","utils","table.html"))
url = "file://"+path
webbrowser.open(url)