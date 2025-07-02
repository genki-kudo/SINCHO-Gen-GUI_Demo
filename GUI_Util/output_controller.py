import streamlit as st
import subprocess
import os
import tempfile
import io
from Bio.PDB import PDBParser
from fractions import Fraction
import py3Dmol
from streamlit.components.v1 import html
import yaml
import shutil
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw



class OutputController:
    def __init__(self):
        # Streamlitのセッション状態を初期化
        if "output_data" not in st.session_state:
            st.session_state.output_data = None



class OutputController:
    def __init__(self):
        st.write("")
    
    def process(self, sub_tab):

        #"Details", ["General", "MD", "SINCHO", "ChemTS", "AAScore"]

        if sub_tab == "General":
            WARN_FLAG = False
            st.title("General Settings")
            st.write("出力ディレクトリを設定し、他タブで結果を可視化できるようにします。")

            st.session_state.output_settings["output_dir"] = st.text_input("出力ディレクトリ", value="6Z0R")
            if not os.path.isdir(os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"])):
                st.warning(f"出力ディレクトリ {st.session_state.output_settings['output_dir']} が存在しません")
                WARN_FLAG = True

            st.session_state.output_settings["MDdir_name"] = st.text_input("MDディレクトリ名", value="01_ConfSamp")
            if not os.path.isdir(os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["MDdir_name"])):
                st.warning(f"ディレクトリ {st.session_state.output_settings['MDdir_name']} が存在しません")
                WARN_FLAG = True
            st.session_state.output_settings["SINCHOdir_name"] = st.text_input("SINCHOディレクトリ名", value="02_MakeDec")
            if not os.path.isdir(os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["SINCHOdir_name"])):
                st.warning(f"ディレクトリ {st.session_state.output_settings['SINCHOdir_name']} が存在しません")
                WARN_FLAG = True
            st.session_state.output_settings["ChemTSdir_name"] = st.text_input("ChemTSディレクトリ名", value="03_CompGen")
            if not os.path.isdir(os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["ChemTSdir_name"])):
                st.warning(f"ディレクトリ {st.session_state.output_settings['ChemTSdir_name']} が存在しません")        
                WARN_FLAG = True          
            st.session_state.output_settings["AAScoredir_name"] = st.text_input("AAScoreディレクトリ名", value="04_DeltaGEst")
            if not os.path.isdir(os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["AAScoredir_name"])):
                st.warning(f"ディレクトリ {st.session_state.output_settings['AAScoredir_name']} が存在しません")
                WARN_FLAG = True
            if not WARN_FLAG:
                st.success("出力ディレクトリの設定が完了しました。")


        elif sub_tab == "MD":
            st.title("MD結果の可視化")
            hr_tabs = st.tabs(["Animation", "Individual"])
            wdir = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["MDdir_name"])
            frames_dir = os.path.join(wdir, "separate_file")
            prot_pdbs = sorted(glob.glob(os.path.join(frames_dir, "prot*.pdb")))
            lig_pdbs = sorted(glob.glob(os.path.join(frames_dir, "lig*.pdb")))
            for line in open(lig_pdbs[0], 'r'):
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    lig_resn = line[17:20].strip()
                    break
            conf_pdbs = sorted(glob.glob(os.path.join(frames_dir, "conf*.pdb")))

            with hr_tabs[0]:
                models = ""
                for i, conf_pdb in enumerate(conf_pdbs):
                    models += f"MODEL {str(i)}\n"
                    models += open(conf_pdb, 'r').read()
                    models += f"ENDMDL\n"
                view = py3Dmol.view(width=700, height=500)
                view.addModelsAsFrames(models)
                view.setStyle({'model':-1}, {"cartoon": {"color": "gray"}})
                view.setStyle({'resn': lig_resn, 'model':-1}, {"stick": {"colorscheme": "greenCarbon"}})
                view.zoomTo({"within": {"distance": 6, "sel": {"resn": lig_resn}}})
                view.animate({'loop': 'forward', 'interval': None})
                html(view._make_html(), height=500, width=700)
                    

            with hr_tabs[1]:
                frame_id = st.slider("フレーム番号を選択", min_value=0, max_value=len(prot_pdbs)-1, value=0)

                view = py3Dmol.view(width=700, height=500)
                view.addModel(open(prot_pdbs[frame_id], 'r').read(), "pdb")
                view.addModel(open(lig_pdbs[frame_id], 'r').read().replace("ATOM  ","HETATM"), "pdb")
                view.setStyle({'model':0}, {"cartoon": {"color": "gray"}})
                view.setStyle({'model':1}, {"stick": {"colorscheme": "greenCarbon"}})
                view.zoomTo({'model':0})
                html(view._make_html(), height=500, width=700)






        
        elif sub_tab == "SINCHO":
            st.title("SINCHO結果の可視化")
            hr_tabs = st.tabs(["Summary", "Individual"])

            sincho_df = pd.DataFrame([], columns=['trajectory_num', 'sincho_rank','anchor atom','predicted mw', 'predicted logp', 'pocket id'])

            trajectories = glob.glob(os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["SINCHOdir_name"], "trajectory_*"))
            index = 0
            for tra in trajectories:
                tnum = os.path.basename(tra).split("_")[1]
                sincho_result_file = os.path.join(tra, "sincho_result.yaml")
                if os.path.isfile(sincho_result_file):
                    with open(sincho_result_file, 'r') as f:
                        sincho_result = yaml.safe_load(f)
                    for result in sincho_result['SINCHO_result']:
                        rnum = result.split("_")[1]
                        anchor_atom = sincho_result['SINCHO_result'][result]["atom_num"].split("_")[-1]
                        predicted_mw = round(float(sincho_result['SINCHO_result'][result]["mw"]),5)
                        predicted_logp = round(float(sincho_result['SINCHO_result'][result]["logp"]),5)
                        pocket = sincho_result['SINCHO_result'][result]["atom_num"].split("_")[0]

                        sincho_df.loc[index] = [tnum, rnum, anchor_atom, predicted_mw, predicted_logp, pocket]
                        index+=1
            sincho_df = sincho_df.sort_values(by=['trajectory_num', 'sincho_rank']).reset_index(drop=True)
            
            with hr_tabs[0]:
                st.write("各トラジェクトリのsincho_result.yamlを読み込み、csv形式で表示します。")
                col1, col2 = st.columns([3,1])
                with col1:
                    st.dataframe(sincho_df, use_container_width=True)
                with col2:
                    lig = os.path.join(trajectories[0], "lig_"+trajectories[0].split("/")[-1].split("_")[-1]+".pdb")
                    self._ligfile_3dview(lig)

            with hr_tabs[1]:
                st.write("個別のSINCHO結果について可視化")
                tra = st.selectbox("トラジェクトリを選択", sincho_df['trajectory_num'].unique())
                ran = st.selectbox("ランクを選択", sincho_df[sincho_df['trajectory_num'] == tra]['sincho_rank'].unique())
                sincho_row = sincho_df[(sincho_df['trajectory_num'] == tra) & (sincho_df['sincho_rank'] == ran)]
                st.write(sincho_row)
                self._sincho_3dview(sincho_row)


        elif sub_tab == "ChemTS":
            st.title("ChemTS結果の可視化")
            hr_tabs = st.tabs(["Summary", "Reward Transition", "Components"])

        elif sub_tab == "AAScore":
            st.title("スコアリング結果の可視化")
            hr_tabs = st.tabs(["General", "Ligand Efficiency"])

            with hr_tabs[0]:
                states = "General"

                csv = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["AAScoredir_name"], "all.csv")
                df_gen = pd.read_csv(csv)
                sdf = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["AAScoredir_name"], "all.sdf")
                suppl = Chem.SDMolSupplier(sdf, removeHs=False)
                mol_indices = list(range(len(suppl)))

                st.session_state.output_settings["mode"] = st.selectbox("表示モードを選択", ["Simple CSV Summary", "2D Summary", "3D View"], key="AAScore_mode")
                if st.session_state.output_settings["mode"] == "Simple CSV Summary":
                    st.write("スコアリング結果をCSV形式で出力します。")

                    st.dataframe(df_gen[["trajectory_num", "rank_num", "lead_num", "AAScore"]], use_container_width=True)
                elif st.session_state.output_settings["mode"] == "2D Summary":
                    st.write("スコアリング結果を2D形式で表示します。")
                    c1, c2 = st.columns([1, 1])
                    with c1:
                        start = st.number_input("Start index", min_value=0, max_value=df_gen.shape[0], value=0, step=1, key="start_index_gen")
                    with c2:
                        end = st.number_input("End index",min_value=start,max_value=min(start + 50, df_gen.shape[0]),value=min(start + 50, df_gen.shape[0]),step=1, key="end_index_gen")
                    if st.button(f"選択範囲: {start} - {end}で２次元描画する", key="gen_2dview_button"):
                        img = self._summary_2Dview(suppl, int(start), int(end), states)

                elif st.session_state.output_settings["mode"] == "3D View":
                    st.write("スコアリング結果を3D形式で表示します。")
                    idx = st.number_input(
                        "表示したい構造のインデックスを入力",
                        min_value=min(mol_indices),
                        max_value=max(mol_indices),
                        value=0,
                        step=1,
                        key="general_3dview_idx"
                    )
                    st.write(df_gen.iloc[idx][["trajectory_num", "rank_num", "lead_num", "AAScore"]])

                    mol = suppl[idx]
                    if mol is None:
                        st.warning(f"{idx} 番目の分子は読み込みできませんでした。")
                    else:
                        trajectory_num = mol.GetProp("trajectory_num")
                        sdir = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["SINCHOdir_name"])
                        prot_pdb = os.path.join(sdir, trajectory_num, f"prot_{trajectory_num.split('_')[-1]}.pdb")
                        self._comp_3dview(mol, prot_pdb)


            with hr_tabs[1]:
                states = "Ligand Efficiency"

                csv = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["AAScoredir_name"], "all_le.csv")
                df_gen = pd.read_csv(csv)
                sdf = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["AAScoredir_name"], "all_le.sdf")
                suppl = Chem.SDMolSupplier(sdf, removeHs=False)
                mol_indices = list(range(len(suppl)))

                st.session_state.output_settings["mode"] = st.selectbox("表示モードを選択", ["Simple CSV Summary", "2D Summary", "3D View"], key="le_mode")
                if st.session_state.output_settings["mode"] == "Simple CSV Summary":
                    st.write("スコアリング結果をCSV形式で出力します。")

                    st.dataframe(df_gen[["trajectory_num", "rank_num", "lead_num", "AAScore_LE"]], use_container_width=True)
                    #st.dataframe(df_gen[["trajectory_num", "rank_num", "lead_num", "AAScore"]], use_container_width=True)
                    
                elif st.session_state.output_settings["mode"] == "2D Summary":
                    st.write("スコアリング結果を2D形式で表示します。")
                    c1, c2 = st.columns([1, 1])
                    with c1:
                        start = st.number_input("Start index", min_value=0, max_value=df_gen.shape[0], value=0, step=1, key="start_index_le")
                    with c2:
                        end = st.number_input("End index",min_value=start,max_value=min(start + 50, df_gen.shape[0]),value=min(start + 50, df_gen.shape[0]),step=1, key="end_index_le")
                    if st.button(f"選択範囲: {start} - {end}で２次元描画する", key="le_2dview_button"):
                        img = self._summary_2Dview(suppl, int(start), int(end), states)

                elif st.session_state.output_settings["mode"] == "3D View":
                    st.write("スコアリング結果を3D形式で表示します。")
                    idx = st.number_input(
                        "表示したい構造のインデックスを入力",
                        min_value=min(mol_indices),
                        max_value=max(mol_indices),
                        value=0,
                        step=1,
                        key="le_3dview_idx"
                    )
                    #st.write(df_gen.iloc[idx][["trajectory_num", "rank_num", "lead_num", "AAScore"]])
                    st.write(df_gen.iloc[idx][["trajectory_num", "rank_num", "lead_num", "AAScore_LE"]])


                    mol = suppl[idx]
                    if mol is None:
                        st.warning(f"{idx} 番目の分子は読み込みできませんでした。")
                    else:
                        trajectory_num = mol.GetProp("trajectory_num")
                        sdir = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["SINCHOdir_name"])
                        prot_pdb = os.path.join(sdir, trajectory_num, f"prot_{trajectory_num.split('_')[-1]}.pdb")
                        self._comp_3dview(mol, prot_pdb)







    def _summary_2Dview(self, suppl, start, end, states):
        mols = list(suppl)[start:end]
        ind = start+1
        for mol in mols:
            AllChem.Compute2DCoords(mol)
            mol.SetProp("rank", str(ind))
            ind+=1
        if states=="General":
            img = Draw.MolsToGridImage(mols,molsPerRow=10, useSVG=False, 
            subImgSize=(200,150),legends=["rank"+mol.GetProp('rank')+": "+str(round(float(mol.GetProp('AAScore')),3))+" kcal/mol" for mol in mols], legendFontSize=40)
        elif states=="Ligand Efficiency":
            img = Draw.MolsToGridImage(mols,molsPerRow=10, useSVG=False, 
            subImgSize=(200,150),legends=["rank"+mol.GetProp('rank')+": "+str(round(float(mol.GetProp('AAScore_LE')),3)) for mol in mols], legendFontSize=40)
            #img = Draw.MolsToGridImage(mols,molsPerRow=10, useSVG=False, 
            #subImgSize=(200,150),legends=["rank"+mol.GetProp('rank')+": "+str(round(float(mol.GetProp('AAScore')),3))+" kcal/mol" for mol in mols], legendFontSize=40)
        st.image(img)
        # ここでバイナリに変換
        from io import BytesIO
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)

        st.download_button(
            label="この描画を保存する",
            data=buffer,
            file_name="SINCHO-Gen_2DView.png",
            mime=os.getcwd()
        )
        return
    
    def _sincho_3dview(self, sincho_row):
        # SINCHO結果の3D表示
        trajectory_num = sincho_row['trajectory_num'].values[0]
        dirname = os.path.join(os.getcwd(), st.session_state.output_settings["output_dir"], st.session_state.output_settings["SINCHOdir_name"], "trajectory_"+str(trajectory_num))
        rank_num = sincho_row['sincho_rank'].values[0]
        anchor_atom = sincho_row['anchor atom'].values[0]
        pocket_id = sincho_row['pocket id'].values[0]
        lig_pdb = os.path.join(dirname, f"lig_{trajectory_num}.pdb")
        prot_pdb = os.path.join(dirname, f"prot_{trajectory_num}.pdb")
        poc_pdb = os.path.join(dirname, 'sincho-output', f"{pocket_id}.pdb")

        view = py3Dmol.view(width=700, height=500)
        with open(lig_pdb, 'r') as f:
            lig_block = f.read()
        view.addModel(lig_block, "pdb")
        with open(prot_pdb, 'r') as f:
            prot_block = f.read()
        view.addModel(prot_block, "pdb")
        with open(poc_pdb, 'r') as f:
            poc_block = f.read()
        view.addModel(poc_block, "pdb")
        view.setStyle({'model': 2}, {})
        view.setStyle({'model': 0}, {"stick": {"colorscheme": "greenCarbon"}})
        view.setStyle({'model': 1}, {"cartoon": {"color": "gray"}})
        #view.setStyle({'model': 2}, {"mesh": {"color": "orange", "scale": 0.3}})
        view.addSurface("VDW", {"color": "orange", "opacity": 1, "mesh":False},{'model': 2})

        view.zoomTo({'model': [0,2]})

        poc_block = poc_block.split("\n")

        n = sum(1 for line in poc_block if line.startswith("ATOM  "))
        px = sum(float(line[30:38]) for line in poc_block if line.startswith("ATOM  "))/n
        py = sum(float(line[38:46]) for line in poc_block if line.startswith("ATOM  "))/n
        pz = sum(float(line[46:54]) for line in poc_block if line.startswith("ATOM  "))/n

        for line in lig_block.split("\n"):
            if (line.startswith("HETATM") or line.startswith("ATOM  ")) and line[12:16].strip() == anchor_atom:
                lx = float(line[30:38])
                ly = float(line[38:46])
                lz = float(line[46:54])

        view.setStyle({'model':0, 'atom': anchor_atom}, {"stick": {"colorscheme": "greenCarbon"}, "sphere": {"color": "orange", "scale":0.3}})
        view.addLine({"start": {"x": lx, "y": ly, "z": lz},
                        "end": {"x": px, "y": py, "z": pz},
                        "color": "orange",
                        "linewidth": 7.5})

        html(view._make_html(), height=500, width=700)

    
        
    def _comp_3dview(self, mol, prot):
        molblock = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=700, height=500)
        view.addModel(molblock, "mol")
        view.setStyle({}, {"stick": {"colorscheme": "greenCarbon"}})

        with open(prot, 'r') as f:
            pdb_block = f.read()
        view.addModel(pdb_block, "pdb")
        view.setStyle({'hetflag': False}, {"cartoon": {"color": "gray"}})

        view.zoomTo({"model": 0})
        html(view._make_html(), height=500, width=800)

    def _ligfile_3dview(self, lig):
        with open(lig, 'r') as f:
            pdb_block = f.read()

        # ATOM → HETATM に置換して誤認識を防ぐ
        pdb_block = pdb_block.replace("ATOM  ", "HETATM")

        view = py3Dmol.view(width=250, height=400)
        view.addModel(pdb_block, "pdb")
        view.setStyle({}, {"stick": {"colorscheme": "greenCarbon"}})

        # 1行ずつラベル付け
        for line in pdb_block.splitlines():
            if line.startswith("HETATM"):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                view.addLabel(name, {
                    "position": {"x": x, "y": y, "z": z},
                    "fontColor": "black",
                    "backgroundColor": "white",
                    "fontSize": 15,
                    "inFront": True
                })

        view.zoomTo()
        html(view._make_html(), height=400, width=250)


    def _pdb_3dview(self, pdbfile, zoomres = None):
        #PDBファイルを3次元表示
        pdb_str = pdbfile
        view = py3Dmol.view(height=500, width=800)
        view.addModel(pdb_str, "pdb")

        # 水分子（HOH, WAT）
        view.setStyle({'resn': 'HOH'}, {"sphere": {"color": "skyblue", "radius": 0.3}})
        view.setStyle({'resn': 'WAT'}, {"sphere": {"color": "skyblue", "radius": 0.3}})
        # リガンドなど（水以外の HETATM）
        view.setStyle({'hetflag': True, 'resn': ['HOH', 'WAT'], 'invert': True},
                      {"stick": {"colorscheme": "greenCarbon"}})
        # タンパク質
        view.setStyle({'hetflag': False}, {"cartoon": {"color": "gray"}})
        if zoomres:
            view.zoomTo({'resi': zoomres.split(" ")[1]})
            view.setStyle({'resi':zoomres.split(" ")[1]}, {"stick":{}})
        else:
            view.zoomTo()
        html(view._make_html(), height=500, width=800)
        return
    

    def _residue_parser(self, pdbfile):
        # PDB ファイルを解析して残基一覧を取得
        parser = PDBParser(QUIET=True)
        # BytesIO オブジェクトをテキスト形式に変換
        pdb_text = io.StringIO(pdbfile.getvalue().decode("utf-8"))
        structure = parser.get_structure("uploaded_structure", pdb_text)
        reslist = [
            f"{residue.get_resname()} {residue.id[1]}"
            for model in structure
            for chain in model
            for residue in chain
        ]
        return reslist

    


