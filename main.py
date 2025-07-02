import subprocess
import os





def main():
    script = os.path.join(os.path.dirname(__file__), "gui_controller.py")
    subprocess.run(["streamlit", "run", script])

if __name__ == "__main__":
    main()
