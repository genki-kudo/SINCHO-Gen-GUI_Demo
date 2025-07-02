import streamlit as st
import os
import tempfile
import io
from Bio.PDB import PDBParser
from fractions import Fraction
import py3Dmol
from streamlit.components.v1 import html

from GUI_Util import initialize, tab_manager
from GUI_Util.input_controller import InputController
from GUI_Util.output_controller import OutputController



#Initialize
initialize.init()

#tab settings
main_tab, sub_tab = tab_manager.init_tabs()

if main_tab == "Input":
    InputController().process(sub_tab)

if main_tab == "Output":
    OutputController().process(sub_tab)
        
