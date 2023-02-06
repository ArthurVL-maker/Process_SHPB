# process_SHPB.py
# ----------------------------------------------------------------------------------------------------------
# Processes raw shpb strain gauge data with dispersion correction or simple time shifting analysis, for
# different confinement types.

# This processing algorithm, 'process_SHPB.py', is also available on ORDA (Van Lerberghe, A., Barr, A. D. (2023)),
# see links below. It was inspired by a Matlab script created by Barr (2023), see link below.

# REQUIRES:
# - dispersion.py & dispersion_factors.py: Implementation of Tyas & Pope (2005)
# 'Full correction of first-mode Pochhammer-Chree dispersion effects in experimental pressure bar signals'.

# INPUTS:
# - raw_file: Path to csv file containing oscilloscope data columns Time, Ch1, Ch2, Ch3, Ch4 (string).
# - sample_data: Length (mm), mass (g), and dry mass (g) data for the sample (cell array)
#               [initial_length, mass, dry_mass].
# - confinement: Specify the specimen confinement 'None'/'Ring'/'Reservoir'.
# - correction: Specify dispersion correction, '1', or time shifting, '0', analysis.

# OUTPUT:
# - A folder titled 'Processed_data' with .csv files containing the processed data.

# NOTES:
# - Confining ring uses output bar Young's modulus.

# REFERENCES:
# - Tyas, A., Pope, D.J., (2005). Full correction of first-mode Pochhammer–Chree dispersion effects in experimental
# pressure bar signals. Measurement science and technology, 16(3), p.642.

# MATLAB SOFTWARE:
# - Barr, A. D. (2016) dispersion.m - A Matlab script for phase angle and amplitude correction of pressure bar signals.
# University of Sheffield.
# Software ORDA link: [https://doi.org/10.15131/shef.data.3996876.v1]

# PYTHON SOFTWARE:
# - Van Lerberghe, A., Barr, A. D. (2023) process_SHPB.py - A Python algorithm for stress wave dispersion
# correction in split-Hopkinson pressure bar experiments. University of Sheffield.
# Software ORDA link: [https://doi.org/10.15131/shef.data.21973325]

# AUTHORS:
# Arthur Van Lerberghe (<avanlerberghe1@sheffield.ac.uk>) & Andrew D. Barr (<a.barr@sheffield.ac.uk>).
# ----------------------------------------------------------------------------------------------------------
# Imported modules:
from pathlib import Path
import pandas as pd
import numpy as np
import time

# Imported function:
from dispersion import dispersion


def process_SHPB(raw_file, sample_data, confinement, correction):
    # ----------------------------------------------------------------
    # VARIABLES
    # ----------------------------------------------------------------
    # Time code:
    start = time.time()

    # Results type:
    if correction == 1:
        result_type = '_dispersion_correction'
    else:
        result_type = '_time_shifting'

    # Raw file path:
    file = Path(raw_file)
    print('-' * 60 + '\n' + f"PROCESSING {file.parts[-1].split('.')[0] + result_type} " '\n' + '-' * 60)
    print(f'Original file path:{file}.')

    # Sample:
    sample_initial_length = sample_data[0]  # mm
    sample_mass = sample_data[1]  # g
    sample_dry_mass = sample_data[2]  # g
    sample_diameter = 25  # mm

    # Hop Bars & Gauges:
    # Incident bar - Steel SS-25:
    in_bar_density = 7666  # Bar density, kg/m**3.
    in_bar_diameter = 25  # Bar diameter, mm.
    in_bar_wave_speed = 5422  # Bar wave speed, m/s.
    in_bar_gauge_channel = 0  # Input bar oscilloscope channel.
    in_bar_gauge_factor = 123  # Input bar gauge factor.
    in_bar_gauge_amp = 1  # Input bar signal amplification.
    in_bar_gauge_voltage = 10  # Input bar signal voltage, V.
    in_bar_gauge_offset = 1000  # Distance from strain gauge to sample face, mm.

    # Transmitter bar - Steel SS-25:
    out_bar_density = 7677  # Bar density, kg/m**3.
    out_bar_diameter = 25  # Bar diameter, mm.
    out_bar_wave_speed = 5311  # Bar wave speed, m/s.
    out_bar_gauge_channel = 1  # Output bar oscilloscope channel.
    out_bar_gauge_factor = 120  # Output bar gauge factor.
    out_bar_gauge_amp = 1  # Output bar signal amplification.
    out_bar_gauge_voltage = 10  # Output bar signal voltage, V.
    out_bar_gauge_offset = 250  # Distance from strain gauge to sample face, mm.

    # Confining type:
    confinement_type = str(confinement)  # Specimen confinement 'None'/'Ring'/'Reservoir'

    # Ring:
    if confinement_type == 'Ring':
        ring_outside_diameter = 35  # Outside diameter, mm.
        ring_inside_diameter = 25  # Inside diameter, mm.
        ring_length = 5  # Length, mm.
        ring_gauge_channel = 2  # Oscilloscope channel.
        ring_gauge_factor = 130  # Gauge factor.
        ring_gauge_amp = 1  # Signal amplification.
        ring_gauge_voltage = 10  # Signal voltage, V.
        ring_youngs_modulus = 206  # Young's modulus, GPa.
        print("Confinement type selected:'Ring'")

    # Reservoir:
    elif confinement_type == 'Reservoir':
        reservoir_fluid_wave_speed = 1482  # Wave speed of water 1482 m/s.
        reservoir_thickness = 9.6  # Thickness of fluid annulus at transducer, mm.
        reservoir_gauge_channel = 2  # Reservoir transducer oscilloscope channel.
        reservoir_gauge_factor = 2.90  # Reservoir transducer calibration, mV/MPa.
        reservoir_gauge_voltage = 10  # Reservoir transducer voltage, V.
        print("Confinement type selected:'Reservoir'")

    # None:
    elif confinement_type == 'None':
        print("Confinement type selected:'None'")

    # No confinement type selected:
    else:
        print("No confinement type selected")

    # ----------------------------------------------------------------
    # RAW DATA
    # ----------------------------------------------------------------
    # CSV file format: Relative time, Channel 1, Channel 2, Channel 3 & Channel 4.
    raw_data = pd.read_csv(raw_file, sep=';', skiprows=9, header=None)  # Read csv file.
    time_base = raw_data.iloc[1:3, 0].values  # First two time values, s.
    in_bar_gauge_signal = raw_data.iloc[1:50000, in_bar_gauge_channel + 1].values  # V.
    out_bar_gauge_signal = raw_data.iloc[1:50000, out_bar_gauge_channel + 1].values  # V.

    # ----------------------------------------------------------------
    # AXIAL PROCESSING
    # ----------------------------------------------------------------
    # Strain gauge signals:
    time_step = time_base[1] - time_base[0]  # Oscilloscope time step, s.
    in_bar_gauge_zero = in_bar_gauge_signal[: 1000].mean()  # Mean input bar "no signal" voltage, V.
    out_bar_gauge_zero = out_bar_gauge_signal[: 1000].mean()  # Mean output bar "no signal" voltage, V.

    # Incident bar strains:
    in_bar_strain = ((in_bar_gauge_signal - in_bar_gauge_zero) * 2) / (in_bar_gauge_factor * in_bar_gauge_voltage * in_bar_gauge_amp)  # Input bar strain.
    in_bar_youngs_modulus = (in_bar_wave_speed ** 2) * (in_bar_density / (10 ** 9))  # Input bar Young's modulus.
    in_bar_stress = in_bar_strain * in_bar_youngs_modulus * 1000

    # Transmitter bar strains:
    out_bar_strain = ((out_bar_gauge_signal - out_bar_gauge_zero) * 2) / (out_bar_gauge_factor * out_bar_gauge_voltage * out_bar_gauge_amp)  # Output bar strain.
    out_bar_youngs_modulus = (out_bar_wave_speed ** 2) * (out_bar_density / (10 ** 9))  # Output bar Young's modulus.
    out_bar_stress = out_bar_strain * out_bar_youngs_modulus * 1000

    # Detect pulses:
    trigger_strain = 0.0001  # Absolute strain indicating start of pulse.
    zero_strain = 0.00001  # Absolute strain for "zero" envelope.

    # Incident pulse:
    incident_trigger = np.where(abs(in_bar_strain) > trigger_strain)[0][0]  # Find when the signal first larger than trigger_strain.
    if in_bar_strain[incident_trigger] < 0:
        in_bar_strain = -in_bar_strain  # If incident wave is negative, invert signal.
    incident_start = np.where(np.array(in_bar_strain[0:incident_trigger]) * np.array(in_bar_strain[1:incident_trigger+1]) < 0)[0][-1]  # Find last change of sign before trigger (start of incident pulse).
    incident_end = np.where((np.array(in_bar_strain[incident_start:-1]) * np.array(in_bar_strain[incident_start + 1:])) < 0)[0][1] + incident_start  # Find the next change of sign after trigger (end of incident pulse).
    incident_length = incident_end - incident_start  # Length of the incident pulse.

    # Reflected pulse:
    reflected_trigger = np.where(abs(in_bar_strain[incident_end:]) > trigger_strain/5)[0][1] + incident_end - 1  # Find when signal next has a value larger than trigger_strain.
    # Play with reflected start function - No2 better than No1:
    reflected_start = np.where(abs(in_bar_strain[incident_end:reflected_trigger]) < zero_strain)[0][-1] + incident_end  # Find the last "zero" before the trigger (start of reflected pulse).
    # reflected_start = np.where((np.array(in_bar_strain[incident_end:reflected_trigger - 1]) * np.array(in_bar_strain[incident_end + 1:reflected_trigger])) < 0)[0][0] + incident_end  # Find the next change of sign after the trigger (end of reflected pulse).
    reflected_end = reflected_start + incident_length

    # Transmitted pulse:
    transmitted_trigger = np.where(abs(out_bar_strain) > trigger_strain)[0][0]  # Find when signal first has a value larger than trigger_strain.
    if out_bar_strain[transmitted_trigger] < 0:
        out_bar_strain = -out_bar_strain  # If transmitted wave is negative invert signal.
    transmitted_start = np.where(out_bar_strain[:transmitted_trigger] < 0)[0][-1]  # Find the last "zero" before the trigger (start of transmitted pulse).
    transmitted_end = np.where((np.array(out_bar_strain[transmitted_start:-1]) * np.array(out_bar_strain[transmitted_start + 1:])) < 0)[0][1] + transmitted_start  # Find the next change of sign after trigger (end of transmitted pulse).

    n = 20000  # Desired length of the FFT input (pulse + zero padding).
    signal_cut_off = max(reflected_end, transmitted_end) + incident_length
    signal_cut_off = max(signal_cut_off, n)

    # Time shifting & Dispersion correction:
    fs = 1 / time_step  # Sampling frequency, Hz
    in_bar_youngs_modulus = (in_bar_wave_speed ** 2) * (in_bar_density / 10 ** 9)  # Input bar Young's modulus
    out_bar_youngs_modulus = (out_bar_wave_speed ** 2) * (out_bar_density / 10 ** 9)  # Output bar Young's modulus

    # Create signal cut off-length stress waves:
    in_bar_incident = np.concatenate((np.zeros(incident_start), np.conj(np.array(in_bar_strain[incident_start:reflected_start+1])), np.zeros(signal_cut_off-reflected_start-1)))
    in_bar_reflected = np.concatenate((np.zeros(reflected_start), np.conj(np.array(in_bar_strain[reflected_start:reflected_end+1])), np.zeros(signal_cut_off-reflected_end-1)))
    out_bar_transmitted = np.concatenate((np.zeros(transmitted_start-200), np.conj(np.array(out_bar_strain[transmitted_start-200:transmitted_end+1])), np.zeros(signal_cut_off-transmitted_end-1)))

    # Apply dispersion correction analysis, see documentation for dispersion.py:
    if correction == 1:
        print('Processing with dispersion correction.')
        [in_bar_incident_strain, in_bar_incident_stress] = dispersion(in_bar_incident, fs, in_bar_diameter/2000, in_bar_wave_speed, in_bar_youngs_modulus, in_bar_gauge_offset/1000)
        [in_bar_reflected_strain, in_bar_reflected_stress] = dispersion(in_bar_reflected, fs, in_bar_diameter/2000, in_bar_wave_speed, in_bar_youngs_modulus, -in_bar_gauge_offset/1000)
        [out_bar_transmitted_strain, out_bar_transmitted_stress] = dispersion(out_bar_transmitted, fs, out_bar_diameter/2000, out_bar_wave_speed, out_bar_youngs_modulus, -out_bar_gauge_offset/1000)

    # Apply simple timeshift analysis:
    elif correction == 0:
        print('Processing with simple timeshift analysis, not dispersion correction.')
        in_bar_shift = round(((in_bar_gauge_offset/1000) / in_bar_wave_speed) / time_step)
        out_bar_shift = round(((out_bar_gauge_offset/1000) / out_bar_wave_speed) / time_step)
        in_bar_incident_strain = np.concatenate((np.array(in_bar_incident[-1-in_bar_shift:]), np.array(in_bar_incident[:-in_bar_shift])))
        in_bar_reflected_strain = np.concatenate((np.array(in_bar_reflected[in_bar_shift:]), np.array(in_bar_reflected[:in_bar_shift])))
        out_bar_transmitted_strain = np.concatenate((np.array(out_bar_transmitted[out_bar_shift:]), np.array(out_bar_transmitted[:out_bar_shift])))
        in_bar_incident_stress = in_bar_incident_strain * in_bar_youngs_modulus * 1000
        in_bar_reflected_stress = in_bar_reflected_strain * in_bar_youngs_modulus * 1000
        out_bar_transmitted_stress = out_bar_transmitted_strain * out_bar_youngs_modulus * 1000

    # Specimen interface stresses and strains:
    trigger = np.where(abs(in_bar_incident_strain) > trigger_strain)[0][0]  # Find the new position of incident pulse.
    go = np.where(in_bar_incident_strain[:trigger-1] * in_bar_incident_strain[1:trigger] < 0)[0][-1]  # Find last change of sign before trigger (start of incident pulse).
    stop = np.where(in_bar_incident_strain[go:-1] * in_bar_incident_strain[go+1:] < 0)[0][1] + go  # Find next change of sign after trigger (end of incident pulse).

    in_bar_incident_strain = in_bar_incident_strain[go:stop+1]
    in_bar_reflected_strain = in_bar_reflected_strain[go:stop+1]
    out_bar_transmitted_strain = out_bar_transmitted_strain[go:stop+1]
    in_bar_incident_stress = in_bar_incident_stress[go:stop+1]
    in_bar_reflected_stress = in_bar_reflected_stress[go:stop+1]
    out_bar_transmitted_stress = out_bar_transmitted_stress[go:stop+1]

    stress_factor = ((in_bar_diameter/2)**2) / ((sample_diameter/2)**2)

    sample_front_stress = stress_factor * (in_bar_incident_stress + in_bar_reflected_stress)  # Stress at incident bar specimen face, MPa.
    sample_back_stress = stress_factor * out_bar_transmitted_stress  # Stress at transmitter bar specimen face, MPa.
    sample_mid_stress = (sample_front_stress + sample_back_stress)/2  # Mean axial specimen stress, MPa.

    # Bar displacements, sample strains:
    in_bar_displacement = np.zeros(stop-go+1)
    out_bar_displacement = np.zeros(stop-go+1)
    sample_strain = np.zeros(stop-go+1)
    in_bar_displacement_alt = np.zeros(stop-go+1)
    sample_strain_alt = np.zeros(stop-go+1)

    for i in range(1, stop-go+1):
        in_bar_displacement[i] = in_bar_displacement[i-1] + ((in_bar_incident_strain[i] - in_bar_reflected_strain[i]) * 1000 * time_step * in_bar_wave_speed)  # Cumulative input bar displacement.
        out_bar_displacement[i] = out_bar_displacement[i-1] + (out_bar_transmitted_strain[i] * 1000 * time_step * out_bar_wave_speed)  # Cumulative output bar displacement.
        sample_strain[i] = (in_bar_displacement[i] - out_bar_displacement[i]) / sample_initial_length  # Sample axial strain.
        in_bar_displacement_alt[i] = in_bar_displacement[i-1] + ((in_bar_incident_strain[i]) * 1000 * time_step * in_bar_wave_speed)  # Cumulative input bar displacement, mm.
        sample_strain_alt[i] = 2 * (in_bar_displacement_alt[i] - out_bar_displacement[i]) / sample_initial_length  # Sample axial strain.

    sample_length = sample_initial_length * (1 - sample_strain)  # Sample length.

    # Sample axial strain rate:
    rel_time = time_step * np.arange(0, stop-go+1)  # Relative time, s.
    sample_strain_rate = np.zeros((2, len(sample_strain)))

    for i in range(0, len(sample_strain)-1):
        sample_strain_rate[0, i] = (rel_time[i] + rel_time[i+1])/2
        sample_strain_rate[1, i] = ((sample_length[i] - sample_length[i+1]) / sample_length[i]) / time_step

    # Strain rate:
    sample_strain_rate_1 = sample_strain_rate[0]
    sample_strain_rate_2 = sample_strain_rate[1]

    # ----------------------------------------------------------------
    # RADIAL STRESSES
    # ----------------------------------------------------------------
    sample_initial_volume = sample_initial_length * np.pi * ((in_bar_diameter/2)**2) * 10**(-3)  # Sample initial volume, cm^3.

    if confinement_type == 'Ring':
        print('Processing third input as confining ring strain.' + '\n')
        ring_gauge_signal = raw_data[ring_gauge_channel + 1].iloc[1:n+1]  # Confining ring signal, V.
        maxval = max(abs(np.transpose(ring_gauge_signal)))
        maxloc = np.where(ring_gauge_signal == maxval)[0][0]
        if ring_gauge_signal[maxloc] < 0:
            ring_gauge_signal = -ring_gauge_signal
        ring_gauge_zero = ring_gauge_signal[:1000].mean()
        sample_radial_strain = (ring_gauge_signal[go:stop+1] - ring_gauge_zero) * 4 / (ring_gauge_factor * ring_gauge_voltage * ring_gauge_amp)
        ring_thick_walled_pipe_factor = (((ring_outside_diameter/2)**2) - ((ring_inside_diameter/2)**2)) / (2*(ring_inside_diameter/2)**2)  # Ratio of internal radial stress on the specimen to circumferential stress in the ring
        sample_radial_stress = (ring_thick_walled_pipe_factor * (ring_youngs_modulus * 1000) * sample_radial_strain) * (ring_length / sample_length)  # Radial stress from the ring, MPa.
        sample_volume = sample_initial_volume * (1-sample_strain)  # Soil sample volume, cm^3.
        sample_density = sample_mass / sample_volume  # Sample density, Mg/m^3.
        sample_dry_density = sample_dry_mass / sample_volume  # Sample dry density, Mg/m^3.

    elif confinement_type == 'Reservoir':
        print('Processing third input as water reservoir pressure.' + '\n')
        reservoir_gauge_signal = raw_data[reservoir_gauge_channel + 1].iloc[1:n+1]  # Reservoir pressure transducer signal, V.
        reservoir_gauge_zero = reservoir_gauge_signal[:1000].mean()
        reservoir_stress = ((reservoir_gauge_signal - reservoir_gauge_zero) * 1000) / reservoir_gauge_factor  # Pressure transducer stress, MPa.
        reservoir_transit = (reservoir_thickness / 1000) / reservoir_fluid_wave_speed  # Time for pulse to travel through reservoir fluid, s.
        reservoir_time_steps = round(reservoir_transit / time_step)  # Timeshift in oscilloscope timesteps.
        sample_radial_stress = np.array(reservoir_stress[reservoir_time_steps:-1], reservoir_stress[:reservoir_time_steps - 1])  # Stress at specimen surface.
        sample_radial_stress = sample_radial_stress[go:stop+1]

    else:
        print('No radial stress/strain measurement selected. There will also be no data for volume density.' + '\n')

    # ----------------------------------------------------------------
    # PROCESS & SAVE DATA
    # ----------------------------------------------------------------
    # File format:
    end = '.csv'

    # Folders to save processed data:
    Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/In_bar_results').mkdir(parents=True, exist_ok=True)
    Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Out_bar_results').mkdir(parents=True, exist_ok=True)
    Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Sample_results').mkdir(parents=True, exist_ok=True)
    Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Time_results').mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------------------------
    # Results:
    # ----------------------------------------------------------------
    # In bar results:
    # ----------------------------------------------------------------
    print('-' * 22 + ' INPUT BAR DATA ' + '-' * 22)

    # Save processed in bar constants, in bar signals and in bar results data:
    in_bar_constants = np.transpose(pd.DataFrame((in_bar_density, in_bar_diameter, in_bar_wave_speed, in_bar_gauge_channel,
                                                  in_bar_gauge_factor, in_bar_gauge_amp, in_bar_gauge_voltage, in_bar_gauge_offset,
                                                  in_bar_gauge_zero, in_bar_youngs_modulus)))
    in_bar_signals = np.transpose(pd.DataFrame((in_bar_gauge_signal, in_bar_strain, in_bar_stress)))
    in_bar_results = np.transpose(pd.DataFrame((in_bar_incident, in_bar_reflected, in_bar_incident_strain, in_bar_incident_stress,
                                                in_bar_reflected_strain, in_bar_reflected_stress, in_bar_displacement, in_bar_displacement_alt), dtype=object))

    # Create designated file for processed in bar data:
    new_filename_in_bar_1, new_filename_in_bar_2, new_filename_in_bar_3 = file.parts[-1].split('.')[0] + result_type + "_in_bar_constants" + end, \
                                                                          file.parts[-1].split('.')[0] + result_type + "_in_bar_signals" + end, \
                                                                          file.parts[-1].split('.')[0] + result_type + "_in_bar_results" + end

    # Create new filepath for processed in bar data:
    filepath_in_bar_1, filepath_in_bar_2, filepath_in_bar_3 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/In_bar_results/' + new_filename_in_bar_1), \
                                                              Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/In_bar_results/' + new_filename_in_bar_2), \
                                                              Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/In_bar_results/' + new_filename_in_bar_3)

    # Print new filename created & processed in bar data filepath:
    print(f"New filename: {new_filename_in_bar_1}"), print(f"New filename: {new_filename_in_bar_2}"), print(f"New filename: {new_filename_in_bar_3}")
    print(f"Filepath: {filepath_in_bar_1}"), print(f"Filepath: {filepath_in_bar_2}"), print(f"Filepath: {filepath_in_bar_3}")

    # Save in bar data in corresponding folder:
    in_bar_constants = np.savetxt(filepath_in_bar_1, in_bar_constants, fmt='%s', delimiter=',', header='Density, Diameter, Wave_speed, Gauge_channel, Gauge_factor, Gauge_amp, Gauge_voltage, Gauge_offset, Gauge_zero, Youngs_mod')
    in_bar_signals = np.savetxt(filepath_in_bar_2, in_bar_signals, fmt='%s', delimiter=',', header='Gauge_signal, Strain, Stress')
    in_bar_results = np.savetxt(filepath_in_bar_3, in_bar_results, fmt='%s', delimiter=',', header='Incident, Reflected, Incident_strain, Incident_stress, Reflected_strain, Reflected_stress, Displacement, Displacement_alt')
    print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

    # ----------------------------------------------------------------
    # Out bar results:
    # ----------------------------------------------------------------
    print('-' * 21 + ' OUTPUT BAR DATA ' + '-' * 22)

    # Save processed out bar constants, out bar signals and out bar results data:
    out_bar_constants = np.transpose(pd.DataFrame((out_bar_density, out_bar_diameter, out_bar_wave_speed, out_bar_gauge_channel,
                                                   out_bar_gauge_factor, out_bar_gauge_amp, out_bar_gauge_voltage, out_bar_gauge_offset,
                                                   out_bar_gauge_zero, out_bar_youngs_modulus)))
    out_bar_signals = np.transpose(pd.DataFrame((out_bar_gauge_signal, out_bar_strain, out_bar_stress)))
    out_bar_results = np.transpose(pd.DataFrame((out_bar_transmitted, out_bar_transmitted_strain,
                                                 out_bar_transmitted_stress, out_bar_displacement)))

    # Create designated file for processed out bar data:
    new_filename_out_bar_1, new_filename_out_bar_2, new_filename_out_bar_3 = file.parts[-1].split('.')[0] + result_type + "_out_bar_constants" + end, \
                                                                             file.parts[-1].split('.')[0] + result_type + "_out_bar_signals" + end, \
                                                                             file.parts[-1].split('.')[0] + result_type + "_out_bar_results" + end

    # Created new filepath for processed out bar data:
    filepath_out_bar_1, filepath_out_bar_2, filepath_out_bar_3 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Out_bar_results/' + new_filename_out_bar_1), \
                                                                 Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Out_bar_results/' + new_filename_out_bar_2), \
                                                                 Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Out_bar_results/' + new_filename_out_bar_3)

    # Print new filename created & processed out bar data filepath:
    print(f"New filename: {new_filename_out_bar_1}"), print(f"New filename: {new_filename_out_bar_2}"), print(f"New filename: {new_filename_out_bar_3}")
    print(f"Filepath: {filepath_out_bar_1}"), print(f"Filepath: {filepath_out_bar_2}"), print(f"Filepath: {filepath_out_bar_3}")

    # Save out bar data in corresponding folder:
    out_bar_constants = np.savetxt(filepath_out_bar_1, out_bar_constants, fmt='%s', delimiter=',', header='Density, Diameter, Wave_speed, Gauge_channel, Gauge_factor, Gauge_amp, Gauge_voltage, Gauge_offset, Gauge_zero, Youngs_mod')
    out_bar_signals = np.savetxt(filepath_out_bar_2, out_bar_signals, fmt='%s', delimiter=',', header='Gauge_signal, Strain, Stress')
    out_bar_results = np.savetxt(filepath_out_bar_3, out_bar_results, fmt='%s', delimiter=',', header='Transmitted, Transmitted_strain, Transmitted_stress, Displacement')
    print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

    # ----------------------------------------------------------------
    # Sample results:
    # ----------------------------------------------------------------
    print('-' * 21 + ' SAMPLE BAR DATA ' + '-' * 22)

    # Save processed sample constants and sample results data:
    sample_constants = np.transpose(pd.DataFrame((sample_initial_length, sample_mass, sample_dry_mass, sample_diameter, sample_initial_volume)))
    sample_results = np.transpose(pd.DataFrame((sample_front_stress, sample_back_stress, sample_mid_stress, sample_strain_alt, sample_strain, sample_length, sample_strain_rate_1, sample_strain_rate_2)))

    # Create designated file for processed sample data:
    new_filename_sample_1, new_filename_sample_2 = file.parts[-1].split('.')[0] + result_type + "_sample_constants" + end, \
                                                   file.parts[-1].split('.')[0] + result_type + "_sample_results" + end

    # Created new filepath for processed sample data:
    filepath_sample_1, filepath_sample_2 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Sample_results/' + new_filename_sample_1), \
                                           Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Sample_results/' + new_filename_sample_2),

    # Print new filename created & processed sample data filepath:
    print(f"New filename: {new_filename_sample_1}"), print(f"New filename: {new_filename_sample_2}")
    print(f"Filepath: {filepath_sample_1}"), print(f"Filepath: {filepath_sample_2}")

    # Save sample data in corresponding folder:
    sample_constants = np.savetxt(filepath_sample_1, sample_constants, fmt='%s', delimiter=',', header='Initial_length, Mass, Dry_mass, Diameter, Initial_volume')
    sample_results = np.savetxt(filepath_sample_2, sample_results, fmt='%s', delimiter=',', header='Front_stress, Back_stress, Mid_stress, Strain_alt, Strain, Length, Strain_rate_1, Strain_rate_2')
    print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

    # ----------------------------------------------------------------
    # Time results:
    # ----------------------------------------------------------------
    print('-' * 24 + ' TIME DATA ' + '-' * 25)

    # Save processed time data:
    time_results = pd.DataFrame(rel_time)

    # Create designated file for processed time data:
    new_filename_time = file.parts[-1].split('.')[0] + result_type + "_time_results" + end

    # Created new filepath for processed time data:
    filepath_time = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Time_results/' + new_filename_time)

    # Print new filename created & processed time data filepath:
    print(f"New filename: {new_filename_time}")
    print(f"Filepath: {filepath_time}")

    # Save time data in corresponding folder:
    time_results = np.savetxt(filepath_time, time_results, fmt='%s', delimiter=',', header='rel_time')
    print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

    if confinement_type == 'Ring':
        # Folder to save processed ring data:
        Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Ring_results').mkdir(parents=True, exist_ok=True)

        # ----------------------------------------------------------------
        # Ring results:
        # ----------------------------------------------------------------
        print('-' * 24 + ' RING DATA ' + '-' * 25)

        # Save processed ring constants and ring signal data:
        ring_constants = np.transpose(pd.DataFrame((ring_outside_diameter, ring_inside_diameter, ring_length, ring_gauge_channel,
                                                    ring_gauge_factor, ring_gauge_amp, ring_gauge_voltage, ring_youngs_modulus,
                                                    ring_gauge_zero, ring_thick_walled_pipe_factor)))
        ring_signal = pd.DataFrame(ring_gauge_signal)

        # Create designated file for processed ring data:
        new_filename_ring_1, new_filename_ring_2 = file.parts[-1].split('.')[0] + result_type + "_ring_constants" + end, \
                                                   file.parts[-1].split('.')[0] + result_type + "_ring_signal" + end

        # Create new filepath for processed ring data:
        filepath_ring_1, filepath_ring_2 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Ring_results/' + new_filename_ring_1), \
                                           Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Ring_results/' + new_filename_ring_2)

        # Print new filename created & processed ring data filepath:
        print(f"New filename: {new_filename_ring_1}"), print(f"New filename: {new_filename_ring_2}")
        print(f"Filepath: {filepath_ring_1}"), print(f"Filepath: {filepath_ring_2}")

        # Save Ring data in corresponding folder:
        ring_constants = np.savetxt(filepath_ring_1, ring_constants, fmt='%s', delimiter=',', header='Outside_diameter, Insider_diameter, Length, Gauge_channel, Gauge_factor, Gauge_amp, Gauge_voltage, Youngs_mod, Gauge_zero, Thick_wall_pipe_factor')
        ring_signal = np.savetxt(filepath_ring_2, ring_signal, fmt='%s', delimiter=',', header='Gauge_signal')
        print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

        # ----------------------------------------------------------------
        # More Sample results:
        # ----------------------------------------------------------------
        print('-' * 21 + ' MORE SAMPLE DATA ' + '-' * 21)

        # Save processed sample properties and sample radial results data:
        sample_properties = np.transpose(pd.DataFrame((sample_volume, sample_density, sample_dry_density)))
        sample_radial_results = np.transpose(pd.DataFrame((sample_radial_strain, sample_radial_stress)))

        # Create designated file for processed more sample data:
        new_filename_sample_3, new_filename_sample_4 = file.parts[-1].split('.')[0] + result_type + "_sample_properties" + end, \
                                                       file.parts[-1].split('.')[0] + result_type + "_sample_radial_results" + end

        # Create new filepath for processed more sample data:
        filepath_sample_3, filepath_sample_4 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Sample_results/' + new_filename_sample_3), \
                                               Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Sample_results/' + new_filename_sample_4)

        # Print new filename created & processed more sample data filepath:
        print(f"New filename: {new_filename_sample_3}"), print(f"New filename: {new_filename_sample_4}")
        print(f"Filepath: {filepath_sample_3}"), print(f"Filepath: {filepath_sample_4}")

        # Save More Sample data in corresponding folder:
        sample_properties = np.savetxt(filepath_sample_3, sample_properties, fmt='%s', delimiter=',', header='Volume, Density, Dry_density')
        sample_radial_results = np.savetxt(filepath_sample_4, sample_radial_results, fmt='%s', delimiter=',', header='Radial_strain, Radial_stress')
        print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

        # Folder to save processed key results data:
        Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Key_results').mkdir(parents=True, exist_ok=True)

    elif confinement_type == 'Reservoir':
        # Folder to save processed reservoir data:
        Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Reservoir_results').mkdir(parents=True, exist_ok=True)

        # ----------------------------------------------------------------
        # Reservoir results:
        # ----------------------------------------------------------------
        print('-' * 22 + ' RESERVOIR DATA ' + '-' * 22)

        # Save processed reservoir constants & reservoir results data:
        reservoir_constants = np.transpose(pd.DataFrame((reservoir_fluid_wave_speed, reservoir_thickness, reservoir_gauge_channel,
                                                         reservoir_gauge_factor, reservoir_gauge_voltage, reservoir_gauge_zero)))

        reservoir_results = np.transpose(pd.DataFrame((reservoir_gauge_signal, reservoir_stress)))

        # Create designated file for processed reservoir data:
        new_filename_reservoir_1, new_filename_reservoir_2 = file.parts[-1].split('.')[0] + result_type + "_reservoir_constants" + end, \
                                                             file.parts[-1].split('.')[0] + result_type + "_reservoir_results" + end

        # Create new filepath for processed reservoir data:
        filepath_reservoir_1, filepath_reservoir_2 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Reservoir_results/' + new_filename_reservoir_1), \
                                                     Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Reservoir_results/' + new_filename_reservoir_2)

        # Print new filename created & processed reservoir data filepath:
        print(f"New filename: {new_filename_reservoir_1}"), print(f"New filename: {new_filename_reservoir_2}")
        print(f"Filepath: {filepath_reservoir_1}"), print(f"Filepath: {filepath_reservoir_2}")

        # Save reservoir data in corresponding folder:
        reservoir_constants = np.savetxt(filepath_reservoir_1, reservoir_constants, fmt='%s', delimiter=',', header='Fluid_wave_speed, Thickness, Gauge_channel, Gauge_factor, Gauge_voltage, Gauge_zero')
        reservoir_results = np.savetxt(filepath_reservoir_2, reservoir_results, fmt='%s', delimiter=',', header='Gauge_signal, Stress')
        print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

        # ----------------------------------------------------------------
        # More Sample results:
        # ----------------------------------------------------------------
        print('-' * 21 + ' MORE SAMPLE DATA ' + '-' * 21)

        # Save processed sample properties and sample radial results data:
        sample_radial_results = np.transpose(pd.DataFrame(sample_radial_stress))

        # Create designated file for processed more sample data:
        new_filename_sample_3 = file.parts[-1].split('.')[0] + result_type + "sample_radial_results" + end

        # Create new filepath for processed more sample data:
        filepath_sample_3 = Path('Processed_data/' + file.parts[-1].split('.')[0] + result_type + '/Sample_results/' + new_filename_sample_3)

        # Print new filename created & processed more sample data filepath:
        print(f"New filename: {new_filename_sample_3}")
        print(f"Filepath: {filepath_sample_3}")

        # Save More Sample data in corresponding folder:
        sample_radial_results = np.savetxt(filepath_sample_3, sample_radial_results, fmt='%s', delimiter=',', header='Radial_stress')
        print('-' * 19 + ' PROCESSING COMPLETED ' + '-' * 19 + '\n')

    # Confirm all processing completed:
    print('-' * 60 + '\n' + f"ALL PROCESSING COMPLETED FOR: {file.parts[-1].split('.')[0] + result_type} " '\n' + '-' * 60)

    # Time code:
    stop = time.time()

    return f'Time to run code: {round(stop-start, 3)}s' + '\n' + '-' * 60
