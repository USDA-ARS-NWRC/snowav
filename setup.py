from setuptools import setup, find_packages
from codecs import open
from os import path

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    # TODO: put package requirements here
]

setup_requirements = [
    # TODO(micahsandusky5): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='snowav',
    version='0.1.0',
    description="Snow and Water Model Analysis and Visualization ",
    long_description=readme + '\n\n' + history,
    author="Mark Robertson",
    author_email='mark.robertson@ars.usda.gov',
    url='https://github.com/roberton-mark/SNOWAV',
    packages=['snowav'
			  ],


    include_package_data=True,
    package_data={'snowav':['./config/CoreConfig.ini']},
    scripts=['./scripts/snow.py'],
    install_requires=requirements,
    license="GPL-3.0",
    zip_safe=False,
    keywords='snowav',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    setup_requires=setup_requirements,

    entry_points={
       'console_scripts': ['snowav = scripts.snow:run',],
    }

)
