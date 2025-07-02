import streamlit as st
import subprocess
import os
import tempfile
import io
from Bio.PDB import PDBParser
from fractions import Fraction
from stmol import showmol
import py3Dmol
from streamlit.components.v1 import html

#Streamlitのセッション状態を初期化する関数
def init_tabs():

    st.sidebar.title("TABS")
    main_tab = st.sidebar.radio("Main Panel",["Input", "Output"])

    # 各タブ内で階層2を定義
    if main_tab == "Input":
        sub_tab = st.sidebar.radio("Details", ["General", "Upload Complex", "Select Hit Ligand", "MD Settings", "SINCHO Settings", "ChemTS Settings", "AAScore Settings", "Summary"])

    elif main_tab == "Output":
        sub_tab = st.sidebar.radio("Details", ["General", "MD", "SINCHO", "ChemTS", "AAScore"])

    return main_tab, sub_tab


