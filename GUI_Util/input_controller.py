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



class InputController:
    def __init__(self):
        st.write("")
    
    def process(self, sub_tab):

        #"General", "Upload Complex", "Select Hit Ligand", "MD Settings", "SINCHO Settings", "ChemTS Settings", "AAScore Settings", "Summary"

        if sub_tab == "General":
            st.title("General Settings")

            st.session_state.general_settings.setdefault("dir_step", "init")

            st.session_state.general_settings["use_num_threads"] = st.number_input(
                "並列スレッド数", value=12, step=1
            )
            st.session_state.general_settings["directory"] = st.text_input(
                "出力ディレクトリ", value="./output"
            )

            wdir = os.path.join(os.getcwd(), st.session_state.general_settings["directory"])

            if st.session_state.general_settings["dir_step"] == "done":
                st.success(f"ディレクトリが設定されました: {wdir}")
                # 次の設定UIなどをここに書ける
                st.write("ここからさらに操作を続けられます")

            elif st.session_state.general_settings["dir_step"] == "init":
                if st.button("上記のディレクトリに設定する"):
                    if os.path.exists(wdir):
                        st.session_state.general_settings["dir_step"] = "confirm"
                        try:
                            st.rerun()
                        except Exception as e:
                            st.experimental_rerun()
                    else:
                        try:
                            os.makedirs(os.path.join(wdir, "99_TMP"), exist_ok=True)
                            st.session_state.general_settings["tmp_dir"] = os.path.join(wdir, "99_TMP")
                            st.session_state.general_settings["dir_step"] = "done"
                            try:
                                st.rerun()
                            except Exception as e:
                                st.experimental_rerun()
                        except Exception as e:
                            st.error(f"作成に失敗しました: {e}")

            elif st.session_state.general_settings["dir_step"] == "confirm":
                st.warning(f"ディレクトリ {wdir} は既に存在します。")
                if st.button("上書きを許可する"):
                    try:
                        os.makedirs(os.path.join(wdir, "99_TMP"), exist_ok=True)
                        st.session_state.general_settings["tmp_dir"] = os.path.join(wdir, "99_TMP")
                        st.session_state.general_settings["dir_step"] = "done"
                        try:
                            st.rerun()
                        except Exception as e:
                            st.experimental_rerun()
                    except Exception as e:
                        st.error(f"作成に失敗しました: {e}")

        if sub_tab == "Upload Complex":
            st.title("複合体立体構造のアップロード")


            uploaded_file = st.file_uploader("複合体立体構造ファイルをアップロードしてください (.pdb)", type=["pdb"])

            if uploaded_file:
                # session_state に記録
                st.session_state.uploaded_pdb_file = uploaded_file
                # 3D可視化
                self._pdb_3dview(uploaded_file)
                
                # 一時フォルダに保存
                tmp_path = os.path.join(st.session_state.general_settings["tmp_dir"], uploaded_file.name)
                with open(tmp_path, "wb") as f:
                    f.write(uploaded_file.getvalue())
                st.session_state.uploaded_pdb_file.path = tmp_path
                
                st.success(f"ファイルがアップロードされました ({tmp_path})")
                st.write("ファイルを変更したい場合は再度アップロードしてください")
                """
                if st.button("アップロードしたファイルをキャンセルする"):
                    st.session_state.uploaded_pdb_file = None
                    uploaded_file = None
                    try:
                        st.rerun()
                    except Exception as e:
                        st.experimental_rerun()
                """
            else:
                if st.session_state.uploaded_pdb_file:
                    # session_state に残っているものを表示
                    self._pdb_3dview(st.session_state.uploaded_pdb_file)
                    st.success(f"既にファイル({st.session_state.uploaded_pdb_file.name})がアップロードされています。")

        if sub_tab == "Select Hit Ligand":
            st.title("ヒット化合物の選択と他残基の確認")

            try:
                # まずPDBがアップロードされているかを確認
                if not st.session_state.uploaded_pdb_file:
                    st.warning("PDB ファイルがアップロードされていません。最初のステップに戻ってください。")
                    st.stop()
                
                st.session_state.residues_list = self._residue_parser(st.session_state.uploaded_pdb_file)
                st.write("アップロードされたファイルの構成を確認してください。")
                
                # session_state で選択値を保持
                if "hit_residue" not in st.session_state or st.session_state.hit_residue not in st.session_state.residues_list:
                    st.session_state.hit_residue = st.session_state.residues_list[0]

                
                selected_residue = st.selectbox(
                    "ヒット化合物となる残基を選択してください(1残基のみ選択可能)",
                    st.session_state.residues_list,
                    index=st.session_state.residues_list.index(st.session_state.hit_residue)
                )
                
                # session_stateに更新
                st.session_state.hit_residue = selected_residue
                
                st.success(f"選択された残基: {st.session_state.hit_residue}")
                self._pdb_3dview(st.session_state.uploaded_pdb_file, zoomres=st.session_state.hit_residue)
                st.success("ヒット化合物が選択されました。次のステップに進んでください。")

            except Exception as e:
                st.error(f"ファイルの解析中にエラーが発生しました: {e}")

        if sub_tab == "MD Settings":
            st.title("MD系構築・計算条件設定")

            if "md_settings" not in st.session_state:
                st.session_state.md_settings = {}
            if "force_field" not in st.session_state.md_settings or len(st.session_state.md_settings["force_field"]) != 4:
                st.session_state.md_settings["force_field"] = ["ff14SB", "gaff2", "TIP3P", "OL3"]
            if "temperature" not in st.session_state.md_settings:
                st.session_state.md_settings["temperature"] = 300
            if "additional_parameters" not in st.session_state.md_settings:
                st.session_state.md_settings["additional_parameters"] = None
            if "box_shape" not in st.session_state.md_settings:
                st.session_state.md_settings["box_shape"] = "rectangular"
            if "box_size" not in st.session_state.md_settings:
                st.session_state.md_settings["box_size"] = 75.0
            if "buffer" not in st.session_state.md_settings:
                st.session_state.md_settings["buffer"] = 15.0
            if "snapshots" not in st.session_state.md_settings:
                st.session_state.md_settings["snapshots"] = None
            if "pr_run_time" not in st.session_state.md_settings:
                st.session_state.md_settings["pr_run_time"] = 10
            if "pr_rec_interval" not in st.session_state.md_settings:
                st.session_state.md_settings["pr_rec_interval"] = 2

            with st.expander("力場パラメータ設定"):
                st.session_state.md_settings["force_field"] = [st.selectbox("タンパク質",["ff14SB", "ff99SB", "ff19SB"], index=["ff14SB", "ff99SB", "ff19SB"].index(st.session_state.md_settings["force_field"][0])),
                                                               st.selectbox("化合物",["gaff2", "gaff"], index=["gaff2", "gaff"].index(st.session_state.md_settings["force_field"][1])),
                                                               st.selectbox("水分子",["TIP3P", "SPC"], index=["TIP3P", "SPC"].index(st.session_state.md_settings["force_field"][2])),
                                                               st.selectbox("RNA (if any)",["OL3", "OL4"], index=["OL3", "OL4"].index(st.session_state.md_settings["force_field"][3]))]
            st.success(f"選択された力場: {st.session_state.md_settings['force_field']}")
            with st.expander("MD系オプション（平衡化過程等は固定値を使用します。今後軽量版平衡化も選択可能にする予定。）"):
                st.session_state.md_settings["temperature"] = st.number_input("温度 (K)(動的変数未実装：現状300K固定です)", value=st.session_state.md_settings["temperature"], step=1)
                st.write("化合物のパラメータファイルの追加アップロード➡無い場合はGasteiger chargeを使用")
                additional_parameters = st.file_uploader("追加のMDパラメータファイル(does not applied)", type=["frcmod", "prep"], accept_multiple_files=True)
                st.session_state.md_settings["additional_parameters"] = additional_parameters
                st.session_state.md_settings["box_shape"] = st.selectbox("ボックス形状", ["rectangular", "cube"], index=["rectangular", "cube"].index(st.session_state.md_settings["box_shape"]))
                if st.session_state.md_settings["box_shape"] == "cube":
                    st.session_state.md_settings["box_size"] = st.number_input("ボックスサイズ (Å)", value=st.session_state.md_settings["box_size"], step=1.0)
                else:
                    st.session_state.md_settings["buffer"] = st.number_input("バッファサイズ (Å)", value=st.session_state.md_settings["buffer"], step=0.1)    

            st.write("Production Runの詳細設定をしてください")
            st.session_state.md_settings["pr_run_time"] = st.number_input("Production Run時間 (ns)", value=st.session_state.md_settings["pr_run_time"], step=1)
            st.session_state.md_settings["pr_rec_interval"] = st.number_input("インターバル (ns)", value= st.session_state.md_settings["pr_rec_interval"], step=1)
            snaps = Fraction(st.session_state.md_settings["pr_run_time"]) / Fraction(st.session_state.md_settings["pr_rec_interval"])
            maximum_steps = st.session_state.md_settings["pr_run_time"]/(0.002/1000000)
            if snaps.denominator==1:
                if int(snaps) > maximum_steps or int(snaps)<1:
                    st.warning(f"1 <= (Production Run時間)/(インターバル) <= {int(maximum_steps)}を満たしてください")
                else:
                    st.success(f"{snaps}コ+1コ(0ns)={snaps+1}コのスナップショットが保存され、以下の処理に使用されます。\n次のステップに進んでください。")
                    st.session_state.md_settings["snapshots"] = int(snaps)
            else:
                st.warning(f"(Production Run時間)/(インターバル)を整数値にしてください")
                st.session_state.md_settings["snapshots"] = None

        if sub_tab == "SINCHO Settings":
            st.title("SINCHOの設定")
            if "p2c_sincho_settings" not in st.session_state:
                st.session_state.p2c_sincho_settings = {}
            if "distance_range" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["distance_range"] = 10.0  # Å
            if "npairs_per_snap" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["npairs_per_snap"] = 10  # ペア数
            if "for_chemts" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["for_chemts"] = 2  # ChemTSに渡す候補数

            if st.session_state.md_settings["snapshots"]:
                st.success(f"{str(st.session_state.md_settings['snapshots']+1)}個のスナップショットを保存して以降の処理を進めます。")
                st.session_state.p2c_sincho_settings["distance_range"] = st.number_input("P2Cのポケット探索範囲(化合物からX[Å]以内)", value=st.session_state.p2c_sincho_settings["distance_range"], step=0.1)
                st.session_state.p2c_sincho_settings["npairs_per_snap"] = st.number_input("SINCHOの予測ペア数(per snap)", value=st.session_state.p2c_sincho_settings["npairs_per_snap"], step=1)
                st.session_state.p2c_sincho_settings["for_chemts"] = st.number_input("各スナップ当たり、上位何候補をChemTSに渡す？", value=st.session_state.p2c_sincho_settings["for_chemts"], step=1)
                if st.session_state.p2c_sincho_settings["npairs_per_snap"] < st.session_state.p2c_sincho_settings["for_chemts"]:
                    st.warning("SINCHOの予測ペア数 > ChemTSに渡すペア数 である必要があります。\n設定を見直してください")
                else:
                    st.success("SINCHOの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")
            else:
                st.warning("スナップショット数の設定が適切ではありません。Production Run設定に戻ってください")

        if sub_tab == "ChemTS Settings":
            st.title("ChemTSv2の設定")
            if "chemts_settings" not in st.session_state:
                st.session_state.chemts_settings = {}
            if "num_chemts_loops" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["num_chemts_loops"] = 4  # 生成の反復回数
            if "c_val" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["c_val"] = 1.0           # C値
            if "threshold_type" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["threshold_type"] = "time"  # 終了条件のタイプ
            if "threshold" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["threshold"] = 0.05      # 終了条件の値 (時間 or 生成数)
            if "function_format" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["function_format"] = "only_sincho"  # 報酬の形式

            st.success(f"{st.session_state.md_settings['snapshots']+1}コのスナップショット×{st.session_state.p2c_sincho_settings['for_chemts']}コのLead Strategyの分だけ、生成を行います。")
            with st.expander("Basic setting"):
                st.session_state.chemts_settings["num_chemts_loops"] = st.number_input("生成の反復回数", value=st.session_state.chemts_settings["num_chemts_loops"], step=1)
                st.session_state.chemts_settings["c_val"] = st.number_input("C値", value=st.session_state.chemts_settings["c_val"], step=0.1)
                st.session_state.chemts_settings["threshold_type"] = st.selectbox("1回の生成の終了条件", options=["time","generation_num"], index=["time","generation_num"].index(st.session_state.chemts_settings["threshold_type"]))
                if st.session_state.chemts_settings["threshold_type"] == "generation_num":
                    st.session_state.chemts_settings["threshold"] = st.number_input("生成数", value=st.session_state.chemts_settings["threshold"], step=1)
                elif st.session_state.chemts_settings["threshold_type"] == "time":
                    st.session_state.chemts_settings["threshold"] = st.number_input("時間 (hour)", value=st.session_state.chemts_settings["threshold"], step=0.01)
                    scaler = st.session_state.chemts_settings["threshold"]*st.session_state.chemts_settings["num_chemts_loops"]*st.session_state.p2c_sincho_settings["npairs_per_snap"]*st.session_state.md_settings["snapshots"]
                    st.write(f"ChemTSのおおよその実行時間: {round(scaler,3)}時間")
                
            with st.expander("Advanced setting"):
                st.write("後々実装予定。それまでは最後の直編集パネルでの設定をお願いします。")
            
            with st.expander("Filter setting"):
                st.write("後々実装予定。それまでは最後の直編集パネルでの設定をお願いします。")
                """
                filter_dict = {"lipinski":[["rule_of_5","rule_of_3"],0],"radical":None,"pubchem":None,"sascore":[3.5],"ring_size":[6],"pains":[["[pains_a]"],0],"donor_acceptor":None}
                for fil,op in filter_dict.items():
                    st.session_state.chemts_settings[f"use_{fil}_filter"] = st.selectbox(f"use_{fil}_filter?", options=["True","False"], index=["True","False"].index(st.session_state.chemts_settings.get(f"use_{fil}_filter")))
                    if st.session_state.chemts_settings[f"use_{fil}_filter"] == "True":
                        if op==None:
                            pass
                        elif len(op)==1:
                            st.session_state.chemts_settings[f"{fil}_threshold"] = st.number_input(f"____{fil}_threshold",value=st.session_state.chemts_settings[f"{fil}_threshold"] ,step=0.1)
                        elif len(op)==2:
                            st.session_state.chemts_settings[f"{fil}_type"] = st.selectbox(f"____{fil}_type", options=op[0], index=op[1])
                """

            with st.expander("Reward setting"):
                st.session_state.chemts_settings["function_format"] = st.selectbox("報酬の形式を選択してください", options=["only_sincho", "cns"], index=["only_sincho", "cns"].index(st.session_state.chemts_settings["function_format"]))
                st.write("各項の形式は固定値を使用します。修正は最後の直編集パネルで行ってください")
                st.write("現状、only_sinchoのみ適用可能。CNSを今後適用予定")
            st.write("ChemTSの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")
            st.success("ChemTSの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")

        if sub_tab == "AAScore Settings":
            st.title("AAScoreの設定")
            if "aascore_settings" not in st.session_state:
                st.session_state.aascore_settings = {}
            if "method" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["method"] = "all"  # スコア計算の方法
            if "num_of_cpd" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["num_of_cpd"] = 50  # ランダム選択する化合物数
            if "reward_cutoff" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["reward_cutoff"] = 1.00  # カットオフのreward値
            if "conf_per_cpd" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["conf_per_cpd"] = 20  # 1化合物当たりのconformation数
            if "max_attempts" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["max_attempts"] = 100  # Embedの最大試行回数
            if "rms_thresh" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["rms_thresh"] = 0.25  # 構造間でpruneするRMS閾値[Å]
            if "protein_range" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["protein_range"] = 13  # 計算時のポケット範囲(ヒットからX[Å]以内の残基)
            if "output_num" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["output_num"] = 5000   # 出力sdfファイルに格納する化合物数



            st.session_state.aascore_settings["method"] = st.selectbox("生成された化合物の内、どれをスコア計算するか？", options=["all", "rand"], index=["all", "rand"].index(st.session_state.aascore_settings["method"]))
            if st.session_state.aascore_settings["method"] == "rand":
                st.session_state.aascore_settings["num_of_cpd"] = st.number_input("ランダムに選択する化合物の数", value=st.session_state.aascore_settings["num_of_cpd"], step=1)
            st.session_state.aascore_settings["reward_cutoff"] = st.number_input("カットオフのreward値", value=st.session_state.aascore_settings["reward_cutoff"], step=0.01, max_value=1.00, min_value=0.00)
            st.session_state.aascore_settings["conf_per_cpd"] = st.number_input("1化合物当たりの最大conformation数", value= st.session_state.aascore_settings["conf_per_cpd"], step=1)
            with st.expander("conformation生成の追加設定"):
                st.session_state.aascore_settings["max_attempts"] = st.number_input("Embedの最大試行回数", value= st.session_state.aascore_settings["max_attempts"], step=1)
                st.session_state.aascore_settings["rms_thresh"] = st.number_input("構造間でpruneするRMS閾値[Å]", value= st.session_state.aascore_settings["rms_thresh"], step=0.01)
            st.session_state.aascore_settings["protein_range"] = st.number_input("計算時のポケット範囲(ヒットからX[Å]以内の残基)", value= st.session_state.aascore_settings["protein_range"], step=1)
            st.session_state.aascore_settings["output_num"] = st.number_input("出力sdfファイルに格納する化合物数", value= st.session_state.aascore_settings["output_num"], step=1)
            st.success("AAScoreの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")

        if sub_tab == "Summary":

            replace_dict = {
                        "__OUTDIR__": str(st.session_state.general_settings["directory"]),
                        "__NUM_THREADS__": str(st.session_state.general_settings["use_num_threads"]),
                        "__INPUT_COMPLEX__": str(os.path.join(st.session_state.general_settings["directory"], "99_TMP", st.session_state.uploaded_pdb_file.name)),
                        "__HIT_RESNAME__": str(st.session_state.hit_residue.split(" ")[0]),
                        "__FORCE_FIELD_PROTEIN__": str(st.session_state.md_settings["force_field"][0]),
                        "__FORCE_FIELD_LIGAND__": str(st.session_state.md_settings["force_field"][1]),
                        "__FORCE_FIELD_WATER__": str(st.session_state.md_settings["force_field"][2]),
                        "__BOX_SHAPE__": str(st.session_state.md_settings["box_shape"]),
                        "__BOX_SIZE__": str(st.session_state.md_settings["box_size"]),
                        "__BUFFER__": str(st.session_state.md_settings["buffer"]),
                        # "__TEMPERATURE__": str(st.session_state.md_settings["temperature"]),
                        "__PRODUCTION_RUNTIME__": str(int(st.session_state.md_settings["pr_run_time"])*1000),
                        "__PRODUCTION_REC_INTERVAL__": str(int(st.session_state.md_settings["pr_rec_interval"])*1000),
                        "__SNAPSHOTS__": str(st.session_state.md_settings["snapshots"]),
                        "__SINCHO_DISTANCE_RANGE__": str(st.session_state.p2c_sincho_settings["distance_range"]),
                        "__SINCHO_NPAIRS_PER_SNAP__": str(st.session_state.p2c_sincho_settings["npairs_per_snap"]),
                        "__SINCHO_FOR_CHEMTS__": str(st.session_state.p2c_sincho_settings["for_chemts"]),
                        "__CHEMTS_NUM_LOOPS__": str(st.session_state.chemts_settings["num_chemts_loops"]),
                        "__CHEMTS_C_VAL__": str(st.session_state.chemts_settings["c_val"]),
                        "__CHEMTS_THRESHOLD_TYPE__": str(st.session_state.chemts_settings["threshold_type"]),
                        "__CHEMTS_THRESHOLD__": str(st.session_state.chemts_settings["threshold"]),
                        #"__CHEMTS_FUNCTION_FORMAT__": str(st.session_state.chemts_settings["function_format"]),
                        "__AASCORE_METHOD__": str(st.session_state.aascore_settings["method"]),
                        "__AASCORE_NUM_OF_CPD__": str(st.session_state.aascore_settings["num_of_cpd"]),
                        "__AASCORE_REWARD_CUTOFF__": str(st.session_state.aascore_settings["reward_cutoff"]),
                        "__AASCORE_CONF_PER_CPD__": str(st.session_state.aascore_settings["conf_per_cpd"]),
                        "__AASCORE_MAX_ATTEMPTS__": str(st.session_state.aascore_settings["max_attempts"]),
                        "__AASCORE_RMS_THRESH__": str(st.session_state.aascore_settings["rms_thresh"]),
                        "__AASCORE_PROTEIN_RANGE__": str(st.session_state.aascore_settings["protein_range"]),
                        "__AASCORE_OUTPUT_NUM__": str(st.session_state.aascore_settings["output_num"]),
                    }


            st.title("設定の確認と実行")
            st.write(f"現在の作業ディレクトリ：{os.getcwd()}に、入力ファイルを作成します。")
            st.session_state.yaml_name = st.text_input("YAMLファイル名", value="conditions_lala.yaml")
            if st.session_state.yaml_content is None:
                with open(os.path.join(os.path.dirname(__file__), "conditions_tmp.yaml"), 'r') as file:
                    yaml_content = file.read()
                    for k, v in replace_dict.items():
                        yaml_content = yaml_content.replace(k, v)
                st.session_state.yaml_content = yaml_content
                for k in replace_dict.keys():
                    if k in st.session_state.yaml_content:
                        st.warning(f"{k} が 指定されていません。他のタブで設定し直してください。")


            st.session_state.yaml_content = st.text_area(
                    "Edit YAML content: 追加で編集したい場合は以下を変更してください。",
                    value=st.session_state.yaml_content,
                    height=600,
                    key="edited_yaml"
                )
            if st.button("yamlの初期化"):
                st.session_state.yaml_content = None
                st.success("YAML内容が初期化されました。")
                try:
                    st.rerun()
                except Exception as e:
                    st.experimental_rerun()
            if st.button("Save YAML to File"):
                overwrite_avoider = 0
                yaml_file = os.path.join(os.getcwd(), st.session_state.yaml_name)
                if not yaml_file.endswith(".yaml"):
                    yaml_file += ".yaml"
                    st.session_state.yaml_name += ".yaml"
                # ファイル名が既に存在する場合は、連番を付けて保存
                while os.path.exists(yaml_file):
                    overwrite_avoider += 1
                    yaml_file = os.path.join(os.getcwd(), f"{st.session_state.yaml_name[:-5]}_{overwrite_avoider}.yaml")


                edited_text = st.session_state.yaml_content
                with open(yaml_file, 'w') as f:
                    f.write(edited_text)
                st.success(f"YAMLファイルを保存しました: {yaml_file}")
                if overwrite_avoider > 0:
                    st.warning(f"ファイル名が重複したため、ファイル名に'_{overwrite_avoider}'を付加しました。")
                try:
                    shutil.copy(yaml_file, os.path.join(st.session_state.general_settings['directory'], st.session_state.general_settings['tmp_dir']))
                    st.success(f"{yaml_file}を{os.path.join(st.session_state.general_settings['directory'], st.session_state.general_settings['tmp_dir'])}にバックアップしました。")
                except:
                    st.error(f"{yaml_file}の{os.path.join(st.session_state.general_settings['directory'], st.session_state.general_settings['tmp_dir'])}に対するバックアップに失敗しました。手動でバックアップしてください。")
                #extend_driver.pyがos.getcwd()になければ、コピーしてくる。
                if not os.path.exists(os.path.join(os.getcwd(), "extend_driver.py")):
                    try:
                        shutil.copy(os.path.join(os.path.dirname(__file__), "extend_driver.py"), os.getcwd())
                        st.success("extend_driver.pyを現在の作業ディレクトリにコピーしました。")
                        st.success("これで全ての設定が完了しました。以下のコマンドを実行して処理を開始してください。")
                        st.code(f"python extend_driver.py {yaml_file.split('/')[-1]}")
                    except Exception as e:
                        st.error(f"extend_driver.pyのコピーに失敗しました: {e}")
                else:
                    st.success("extend_driver.pyは既に現在の作業ディレクトリに存在します。")
                    st.success("これで全ての設定が完了しました。以下のコマンドを実行して処理を開始してください。")
                    st.code(f"python extend_driver.py {yaml_file.split('/')[-1]}")


        








    def _pdb_3dview(self, pdbfile, zoomres = None):
        #PDBファイルを3次元表示
        pdb_str = st.session_state.uploaded_pdb_file.getvalue().decode("utf-8")
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

    


