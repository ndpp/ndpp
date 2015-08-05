try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Nuclear Data Pre-Processor [NDPP]',
    'author': 'Adam Nelson',
    'url': 'URL to get it at.',
    'download_url': 'ndpp.github.io',
    'author_email': 'nelsonag at umich.edu',
    'version': '0.1',
    'install_requires': ['nose, numpy, scipy'],
    'packages': ['pyndpp'],
    'scripts': [],
    'name': 'NDPP'
}

setup(**config)
