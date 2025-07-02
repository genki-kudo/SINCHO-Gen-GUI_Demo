import subprocess
import os





def main():
    script = os.path.join(os.path.dirname(__file__), "gui_controller.py")
    subprocess.run(["streamlit", "run", script, "--server.address", "0.0.0.0"]])

if __name__ == "__main__":
    main()
