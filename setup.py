from setuptools import setup

setup(name='dnacycpv2',
      packages=['dnacycpv2'],
      version='0.0.1dev1',
      python_requires='>3.9.0',
      install_requires=[
      # 'numpy==1.21.5',
      'numpy==1.26.1',
      # 'pandas==1.3.5',
      'pandas==2.1.2',
      # 'tensorflow==2.7.0',
      'tensorflow==2.14.0',
      # 'keras==2.7.0',
      'keras==2.14.0',
      # 'bio==1.3.3',
      'bio==1.7.1',
      'docopt==0.6.2'
      ],
      entry_points={
            'console_scripts': ['dnacycpv2-cli=dnacycpv2.cli:main']
      }
      )
