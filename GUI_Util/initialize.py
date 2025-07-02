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
def init():
    # Streamlit のセッション状態を初期化
    if "uploaded_pdb_file" not in st.session_state:
        st.session_state.uploaded_pdb_file = None
    if "residues_list" not in st.session_state:
        st.session_state.residues_list = []
    if "hit_residue" not in st.session_state:
        st.session_state.hit_residue = None
    if "md_settings" not in st.session_state:
        st.session_state.md_settings = {}
    if "p2c_sincho_settings" not in st.session_state:
        st.session_state.p2c_sincho_settings={}
    if "chemts_settings" not in st.session_state:
        st.session_state.chemts_settings = {}
    if "aascore_settings" not in st.session_state:
        st.session_state.aascore_settings = {}
    if "general_settings" not in st.session_state:
        st.session_state.general_settings = {}
    if "output_settings" not in st.session_state:
        st.session_state.output_settings = {}
    st.session_state.general_settings.setdefault("dir_checker", False)
    st.session_state.general_settings.setdefault("dir_overwrite_confirm", False)

    st.session_state.input_state = False  # 初期化フラグ

    return st
