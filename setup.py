from setuptools import setup

setup(name='dnacycp2',
      packages=['dnacycp2'],
      version='0.0.1dev1',
      python_requires='>3.9.0,<3.12',
      install_requires=[
      'numpy==1.26.1',
      'pandas==2.1.2',
      'tensorflow==2.14.0',
      'keras==2.14.0',
      'bio==1.7.1',
      'docopt==0.6.2'
      ],
      entry_points={
            'console_scripts': ['dnacycp2-cli=dnacycp2.cli:main']
      }
      )
