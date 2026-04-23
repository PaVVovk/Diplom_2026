from input_script import create_bolsig_input_data
import subprocess
import os

temp_name = "temporary_input_data.txt"
bolsig_path = "bolsigminus.exe"

def collect_filename():
    default = "output.dat"
    filename = input("Enter output data filename: ").strip()
    return filename + ".dat" if filename else default

filename = collect_filename()

create_bolsig_input_data(temp_name, filename)

process = subprocess.Popen(
    [bolsig_path],
    stdin=subprocess.PIPE,   # для отправки данных
    stdout=subprocess.PIPE,  # для чтения вывода
    stderr=subprocess.PIPE,
    text=True,               # работаем со строками, а не байтами
    bufsize=1                # построчная буферизация
)

process.stdin.write(temp_name + '\n')
process.stdin.flush()  # отправляем данные
process.wait()

try:
    os.remove(temp_name)
except OSError as e:
    print(f"Temporary files deletion error: {e}")

print(f"BOLSIG+ data file {filename} successfully created!")

