from setuptools import setup
setup(
  name = 'pypolybuilder',
  packages = ['pypolybuilder'], # this must be the same as the name above
  entry_points = {
    "console_scripts": ['pypolybuilder = pypolybuilder.pypolybuilder:main']
  },
  version = '0.1',
  description = 'pypolybuilder',
  author = 'Vitor Horta',
  author_email = 'vitor.horta@gmail.com',
  # url = 'https://github.com/vitorhorta/expertSystemPython', # use the URL to the github repo
  # download_url = 'https://github.com/vitorhorta/es/archive/0.7.tar.gz', # I'll explain this in a second
  keywords = ['pypolybuilder'], # arbitrary keywords
  classifiers = [],
)