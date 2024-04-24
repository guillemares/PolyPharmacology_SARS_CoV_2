import setuptools

with open("README.md", "r", encoding='utf-8') as rna:
    readme = rna.read()

setuptools.setup(
        name = 'aRNAlysis',
        version = '0.1',
        author = 'Guillem Arasa',
        author_email = 'guillemares4@gmail.com',
        description = 'Python packages for 3D structures RNA analysis',
        url = 'https://github.com/guillemares/3D_RNAlysis',
        packages = find_packages(),
        install_requires = [
            'Biotite (0.39.0)',
            'USAlign (Version 20231222)'
            ]
        )


