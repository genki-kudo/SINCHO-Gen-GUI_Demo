import subprocess
import os
import shutil





def main():
    conf_dir = os.path.join(os.path.dirname(__file__), ".streamlit")
    dst_dir = os.path.join(os.path.join(os.getcwd(), ".streamlit"))
    if not os.path.exists(dst_dir):
        shutil.copytree(conf_dir, dst_dir)
    script = os.path.join(os.path.dirname(__file__), "gui_controller.py")
    subprocess.run(["streamlit", "run", script, "--server.address", "0.0.0.0"])

if __name__ == "__main__":
    main()
