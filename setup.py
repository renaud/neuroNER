
from setuptools import setup
setup(
   name = "neuroNER",
   version = "0.1",
   packages = ['neuroner'],
   package_data = {
       'neuroner': ['resources/bluima/neuroner/*'],
   },
)