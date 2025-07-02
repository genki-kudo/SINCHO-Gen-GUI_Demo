from setuptools import setup

setup(name="sincho-gen-gui",
      version="0.1.0",
      py_modules=["gui_controller"],
      entry_points={
          "console_scripts":[
              "sincho-gen-gui = gui_controller:main"
          ]
      },
      install_requires=[
          "streamlit",
            "biopython",
            "py3Dmol",
            "pandas",
            "numpy",
            "rdkit",
      ],)
