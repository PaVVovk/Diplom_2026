import numpy as np
from itertools import product

def get_float(prompt, default):
    val = input(f"{prompt} [default {default}]: ").strip()
    return float(val) if val else default

def get_int(prompt, default):
    val = input(f"{prompt} [default {default}]: ").strip()
    return int(val) if val else default

def get_var(prompt, default):
    val = input(f"{prompt} [default {default}, VAR if variable]: ").strip()
    if val.upper() == "VAR":
        start, end, points = map(float, input("Enter start, end and number of variable points: ").split())
        points = int(points)
        return "VAR", np.linspace(start, end, points)
    if val:
        return float(val), []
    else:
        return default, []

def combination_generator(*arrays):
    """Генератор комбинаций, не хранит всё в памяти"""
    for combo in product(*arrays):
        yield combo

def form_run_str(var_params):
    run_str = ""
    arrays = [x[1] for x in var_params.values()]
    for combo in combination_generator(*arrays):
        run_str += ' '.join(map(str, combo)) + "\n"
    return run_str

def create_bolsig_input_data(data_name, output_name, coll_data_path="LXCat-June2013.txt"):
    print("Ввод параметров CONDITIONS:")

    var_params = {
        'e_n': get_var("Electric field / N (Td)", 1.0),
    }
    
    params = {
        'e_n': var_params['e_n'][0],
        'ang_freq': get_float("Angular field frequency / N (m3/s)", 0.0),
        'cos_eb': get_float("Cosine of E-B field angle", 0.0),
        'gas_temp': get_float("Gas temperature (K)", 300.0),
        'exc_temp': get_float("Excitation temperature (K)", 300.0),
        'trans_energy': get_float("Transition energy (eV)", 0.0),
        'ionization_degree': get_float("Ionization degree", 0.0),
        'plasma_density': get_float("Plasma density (1/m3)", 1e18),
        'ion_charge': get_int("Ion charge parameter", 1),
        'ion_mass_ratio': get_float("Ion/neutral mass ratio", 1.0),
        'ee_momentum': get_int("e-e momentum effects (0=No;1=Yes)", 1),
        'energy_sharing': get_int("Energy sharing (1=Equal;2=One takes all)", 1),
        'growth': get_int("Growth (1=Temporal;2=Spatial;3=Not included;4=Grad-n)", 1),
        'maxwell_energy': get_float("Maxwellian mean energy (eV)", 0.0),
        'grid_points': get_int("# of grid points", 100),
        'manual_grid': get_int("Manual grid (0=No;1=Linear;2=Parabolic)", 0),
        'max_energy': get_float("Manual maximum energy (eV)", 200.0),
        'precision': get_float("Precision", 1e-10),
        'convergence': get_float("Convergence", 1e-4),
        'max_iter': get_int("Maximum # of iterations", 2000),
    }

    run_str = form_run_str(var_params)
    
    output = f"""READCOLLISIONS
{coll_data_path} / File
Ar / Species
1 / Extrapolate: 0= No 1= Yes

CONDITIONS
{params['e_n']} / Electric field / N (Td)
{params['ang_freq']} / Angular field frequency / N (m3/s)
{params['cos_eb']} / Cosine of E-B field angle
{params['gas_temp']} / Gas temperature (K)
{params['exc_temp']} / Excitation temperature (K)
{params['trans_energy']} / Transition energy (eV)
{params['ionization_degree']} / Ionization degree
{params['plasma_density']} / Plasma density (1/m3)
{params['ion_charge']} / Ion charge parameter
{params['ion_mass_ratio']} / Ion/neutral mass ratio
{params['ee_momentum']} / e-e momentum effects: 0=No; 1=Yes*
{params['energy_sharing']} / Energy sharing: 1=Equal*; 2=One takes all
{params['growth']} / Growth: 1=Temporal*; 2=Spatial; 3=Not included; 4=Grad-n expansion
{params['maxwell_energy']} / Maxwellian mean energy (eV)
{params['grid_points']} / # of grid points
{params['manual_grid']} / Manual grid: 0=No; 1=Linear; 2=Parabolic
{params['max_energy']} / Manual maximum energy (eV)
{params['precision']} / Precision
{params['convergence']} / Convergence
{params['max_iter']} / Maximum # of iterations

RUN
{run_str}

SAVERESULTS
{output_name} / File
2 / Format: 1=Run by run; 2=Combined; 3=E/N; 4=Energy; 5=SIGLO; 6=PLASIMO
0 / Conditions: 0=No; 1=Yes
1 / Transport coefficients: 0=No; 1=Yes
0 / Rate coefficients: 0=No; 1=Yes
0 / Reverse rate coefficients: 0=No; 1=Yes
0 / Energy loss coefficients: 0=No; 1=Yes
0 / Distribution function: 0=No; 1=Yes
0 / Skip failed runs: 0=No; 1=Yes"""
    
    filename = data_name
    with open(filename, 'w') as f:
        f.write(output)
    #print(f"Файл {filename} создан")
